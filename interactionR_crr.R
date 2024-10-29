library(msm)

interactionR_crr <- function(model, exposure_pos = c(), ci.type = "delta", ci.level = 0.95, 
                             em = TRUE, recode = FALSE) {
  # Check model type
  if (!inherits(model, c("crr"))) {
    stop("Model type not supported. Please provide crr model.")
  }
  
  # Check exposure names
  check_arguments <- function(model, exposure_pos) {
    if (length(exposure_pos) < 2) {
      stop("Please provide exactly two exposure position.")
    }
  }
  check_arguments(model, exposure_pos)
  
  alpha <- 1 - ci.level
  z <- qnorm(1 - alpha / 2)
  
  # Extract coefficients based on model type
  if (inherits(model, "crr")) {
    coef_model <- model$coef
    var_model <- model$var
    
    coef_model <- fg$coef
    var_model <- fg$var
    
  } else {
    coef_model <- coef(model)
    var_model <- vcov(model)
  }
  
  beta1 <- rownames(summary(fg)$coef)[exposure_pos[1]]
  beta2 <- rownames(summary(fg)$coef)[exposure_pos[2]]
  beta3 <- rownames(summary(fg)$coef)[exposure_pos[3]]
  
  beta1_1 <- exposure_pos[1]
  beta2_1 <- exposure_pos[2]
  beta3_1 <- exposure_pos[3]
  
  varNames <- c(beta1, beta2, beta3)
  
  b1 <- as.numeric(coef_model[beta1])
  b2 <- as.numeric(coef_model[beta2])
  b3 <- as.numeric(coef_model[beta3])
  se_vec <- summary(fg)$coef[,3]
  v1 <- se_vec[beta1]^2
  v2 <- se_vec[beta2]^2
  v3 <- se_vec[beta3]^2
  pvals <- summary(fg)$coef[,"p-value"]
  v_cov1 <- var_model[exposure_pos, exposure_pos, drop = FALSE]
  cov12 <- var_model[beta1_1, beta2_1]
  cov13 <- var_model[beta1_1, beta3_1]
  cov23 <- var_model[beta2_1, beta3_1]
  v123 <- v1 + v2 + v3 + (2 * (cov12 + cov13 + cov23))
  v12 <- v1 + v2 + (2 * (cov12))
  v13 <- v1 + v3 + (2 * (cov13))
  v23 <- v2 + v3 + (2 * (cov23))
  
  summaryM = summary(model)
  
  OR00 <- 1
  OR10 <- summaryM$coef[beta1, 2]
  l1 <- summaryM$conf.int[beta1, 3]
  u1 <- summaryM$conf.int[beta1, 4]
  p.OR10 <- summaryM$coef[beta1, 5]
  
  OR01 <- summaryM$coef[beta2, 2]
  l2 <- summaryM$conf.int[beta2, 3]
  u2 <- summaryM$conf.int[beta2, 4]
  p.OR01 <-  summaryM$coef[beta2, 5]
  
  OR11 <- as.numeric(exp(b1 + b2 + b3))
  l3 <- exp(b1 + b2 + b3 - z * sqrt(v123))
  u3 <- exp(b1 + b2 + b3 + z * sqrt(v123))
  q1 <- abs(log(OR11)/sqrt(v123))
  p.OR11 <- exp(-0.717 * q1 - 0.416 * q1^2)
  
  OR_X1 <- as.numeric(exp(b2 + b3))
  CI.ll_OR_X1 <- exp(b2 + b3 - z * sqrt(v23))
  CI.ul_OR_X1 <- exp(b2 + b3 + z * sqrt(v23))
  q2 <- abs(log(OR_X1)/sqrt(v23))
  p.OR_X1 <- exp(-0.717 * q2 - 0.416 * q2^2)
  
  OR_A1 <- as.numeric(exp(b1 + b3))
  CI.ll_OR_A1 <- exp(b1 + b3 - z * sqrt(v13))
  CI.ul_OR_A1 <- exp(b1 + b3 + z * sqrt(v13))
  q3 <- abs(log(OR_A1)/sqrt(v13))
  p.OR_A1 <- exp(-0.717 * q3 - 0.416 * q3^2)
  
  OR_M <- summaryM$coef[beta3, 2]
  CI.ll_OR_M <- summaryM$conf.int[beta3, 3]
  CI.ul_OR_M <- summaryM$conf.int[beta3, 4]
  p.OR_M <- summaryM$coef[beta3, 5]
  if (ci.type == "mover") {
    RERI <- OR11 - OR01 - OR10 + 1
    r12 <- (v1 + cov12 + cov13)/sqrt(v1 * v123)
    r13 <- (cov12 + v2 + cov23)/sqrt(v2 * v123)
    r23 <- cov12/sqrt(v1 * v2)
    p1 <- (OR11 - l3)^2 + (u1 - OR10)^2 + (u2 - OR01)^2
    p2 <- 2 * r12 * (OR11 - l3) * (u1 - OR10)
    p3 <- 2 * r13 * (OR11 - l3) * (u2 - OR01)
    p4 <- 2 * r23 * (u1 - OR10) * (u2 - OR01)
    p5 <- p1 - p2 - p3 + p4
    p6 <- p5^0.5
    L <- 1 + OR11 - OR10 - OR01 - p6
    k1 <- (u3 - OR11)^2 + (OR10 - l1)^2 + (OR01 - l2)^2
    k2 <- 2 * r12 * (u3 - OR11) * (OR10 - l1)
    k3 <- 2 * r13 * (u3 - OR11) * (OR01 - l2)
    k4 <- 2 * r23 * (OR10 - l1) * (OR01 - l2)
    k5 <- (k1 - k2 - k3 + k4)^0.5
    U <- 1 + OR11 - OR10 - OR01 + k5
    p.RERI <- NA
    theta1 <- 1/exp(b1 + b2 + b3)
    theta2 <- 1/exp(b2 + b3)
    theta3 <- 1/exp(b1 + b3)
    AP <- theta1 - theta2 - theta3 + 1
    APr12 <- (cov12 + cov13 + v2 + (2 * cov23) + v3)/sqrt(v23 * 
                                                            v123)
    APr13 <- (v1 + cov12 + (2 * cov13) + cov23 + v3)/sqrt(v13 * 
                                                            v123)
    APr23 <- (cov12 + cov23 + cov13 + v3)/sqrt(v23 * v13)
    APl1 <- theta1 * exp(-z * sqrt(v123))
    APu1 <- theta1 * exp(z * sqrt(v123))
    APl2 <- theta2 * exp(-z * sqrt(v23))
    APu2 <- theta2 * exp(z * sqrt(v23))
    APl3 <- theta3 * exp(-z * sqrt(v13))
    APu3 <- theta3 * exp(z * sqrt(v13))
    APp1 <- (theta1 - APl1)^2 + (APu2 - theta2)^2 + (APu3 - 
                                                       theta3)^2
    APp2 <- 2 * APr12 * (theta1 - APl1) * (APu2 - theta2)
    APp3 <- 2 * APr13 * (theta1 - APl1) * (APu3 - theta3)
    APp4 <- 2 * APr23 * (APu2 - theta2) * (APu3 - theta3)
    APp5 <- APp1 - APp2 - APp3 + APp4
    APp6 <- APp5^0.5
    APL <- 1 + theta1 - theta2 - theta3 - APp6
    APk1 <- (APu1 - theta1)^2 + (theta2 - APl2)^2 + (theta3 - 
                                                       APl3)^2
    APk2 <- 2 * APr12 * (APu1 - theta1) * (theta2 - APl2)
    APk3 <- 2 * APr13 * (APu1 - theta1) * (theta3 - APl3)
    APk4 <- 2 * APr23 * (theta2 - APl2) * (theta3 - APl3)
    APk5 <- (APk1 - APk2 - APk3 + APk4)^0.5
    APU <- 1 + theta1 - theta2 - theta3 + APk5
    p.AP <- NA
    SItheta1 <- log((exp(b1 + b2 + b3) - 1))
    SItheta2 <- log((exp(b1) + exp(b2) - 2))
    lnSI <- SItheta1 - SItheta2
    SI <- exp(lnSI)
    vSItheta1 <- (exp(b1 + b2 + b3)/(exp(b1 + b2 + b3) - 
                                       1))^2 * v123
    vSItheta2 <- ((exp(2 * b1) * v1) + (exp(2 * b2) * v2) + 
                    (2 * exp(b1 + b2) * cov12))/(exp(b1) + exp(b2) - 
                                                   2)^2
    SIl1 <- SItheta1 - z * sqrt(vSItheta1)
    SIu1 <- SItheta1 + z * sqrt(vSItheta1)
    SIl2 <- SItheta2 - z * sqrt(vSItheta2)
    SIu2 <- SItheta2 + z * sqrt(vSItheta2)
    SIr <- ((exp(b1) * (v1 + cov12 + cov13)) + (exp(b2) * 
                                                  (cov12 + v2 + cov23)))/sqrt(v123 * ((exp(2 * b1) * 
                                                                                         v1) + (exp(2 * b2) * v2) + (2 * exp(b1 + b2) * cov12)))
    lnSIL <- (SItheta1 + (-SItheta2)) - sqrt((SItheta1 - 
                                                SIl1)^2 + ((-SItheta2) - (-SIl2))^2 + (2 * SIr * 
                                                                                         (SItheta1 - SIl1) * ((-SItheta2) - (-SIl2))))
    lnSIU <- (SItheta1 + (-SItheta2)) + sqrt((SIu1 - SItheta1)^2 + 
                                               ((-SIu2) - (-SItheta2))^2 + (2 * SIr * (SIu1 - SItheta1) * 
                                                                              ((-SIu2) - (-SItheta2))))
    SIL <- exp(lnSIL)
    SIU <- exp(lnSIU)
    p.SI <- NA
  }
  else if (ci.type == "delta") {
    RERI <- OR11 - OR01 - OR10 + 1
    se_RERI <- deltamethod(g = ~exp(x1 + x2 + x3) - exp(x1) - 
                             exp(x2) + 1, mean = c(b1, b2, b3), cov = v_cov1)
    L <- RERI - z * se_RERI
    U <- RERI + z * se_RERI
    p.RERI <- 1 - pnorm(RERI/se_RERI)
    AP <- RERI/OR11
    se_AP <- deltamethod(g = ~(exp(x1 + x2 + x3) - exp(x1) - 
                                 exp(x2) + 1)/exp(x1 + x2 + x3), mean = c(b1, b2, 
                                                                          b3), cov = v_cov1)
    APL <- AP - z * se_AP
    APU <- AP + z * se_AP
    p.AP <- 1 - pnorm(abs(AP)/se_AP)
    lnSI <- log((exp(b1 + b2 + b3) - 1)) - log((exp(b1) + 
                                                  exp(b2) - 2))
    SI <- exp(lnSI)
    se_SI <- deltamethod(g = ~log((exp(x1 + x2 + x3) - 1)) - 
                           log((exp(x1) + exp(x2) - 2)), mean = c(b1, b2, b3), 
                         cov = v_cov1)
    SIL <- exp(lnSI - z * se_SI)
    SIU <- exp(lnSI + z * se_SI)
    p.SI <- 1 - plnorm(exp(lnSI/se_SI))
  }
  else {
    stop("Argument 'ci.type' must be 'delta' or 'mover' ")
  }
  d <- data.frame(Measures = c("OR00", "OR01", "OR10", "OR11", 
                               paste("OR(", beta2, " on outcome [", beta1, "==0]", 
                                     sep = ""), paste("OR(", beta2, " on outcome [", 
                                                      beta1, "==1]", sep = ""), "Multiplicative scale", 
                               "RERI"), Estimates = c(OR00, OR01, OR10, OR11, OR01, 
                                                      OR_X1, OR_M, RERI), CI.ll = c(NA, l2, l1, l3, l2, CI.ll_OR_X1, 
                                                                                    CI.ll_OR_M, L), CI.ul = c(NA, u2, u1, u3, u2, CI.ul_OR_X1, 
                                                                                                              CI.ul_OR_M, U), p = c(NA, p.OR01, p.OR10, p.OR11, p.OR01, 
                                                                                                                                    p.OR_X1, p.OR_M, p.RERI))
  rownames(d) <- NULL
  if (!em) {
    d <- data.frame(Measures = c("OR00", "OR01", "OR10", 
                                 "OR11", paste("OR(", beta2, " on outcome [", beta1, 
                                               "==0]", sep = ""), paste("OR(", beta2, " on outcome [", 
                                                                        beta1, "==1]", sep = ""), paste("OR(", beta1, 
                                                                                                        " on outcome [", beta2, "==0]", sep = ""), paste("OR(", 
                                                                                                                                                         beta1, " on outcome [", beta2, "==1]", sep = ""), 
                                 "Multiplicative scale", "RERI", "AP", "SI"), Estimates = c(OR00, 
                                                                                            OR01, OR10, OR11, OR01, OR_X1, OR10, OR_A1, OR_M, 
                                                                                            RERI, AP, SI), CI.ll = c(NA, l2, l1, l3, l2, CI.ll_OR_X1, 
                                                                                                                     l1, CI.ll_OR_A1, CI.ll_OR_M, L, APL, SIL), CI.ul = c(NA, 
                                                                                                                                                                          u2, u1, u3, u2, CI.ul_OR_X1, u1, CI.ul_OR_A1, CI.ul_OR_M, 
                                                                                                                                                                          U, APU, SIU), p = c(NA, p.OR01, p.OR10, p.OR11, 
                                                                                                                                                                                              p.OR01, p.OR_X1, p.OR10, p.OR_A1, p.OR_M, p.RERI, 
                                                                                                                                                                                              p.AP, p.SI))
    rownames(d) <- NULL
  }
  
  d$p = round(d$p, 5)
  ir <- list(dframe = d, exp_names = c(beta1, beta2), analysis = em, 
             call = model$call)
  attr(ir, "class") <- "interactionR"
  invisible(ir)
  
}
