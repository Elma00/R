
# interactionR for Competitive Risk Models

This repository contains functions and demo scripts for extending the **interactionR** package to support **competitive risk models**, particularly using the `crr()` function from the **cmprsk** package. These tools allow for modeling interactions within cumulative incidence functions, providing more flexibility and depth in competitive risk analysis.

## Features

- **interactionR_crr.R**: A function that enhances the interactionR functionality to work with `crr()` from the cmprsk package.
- **interactionR_crr_demo.R**: A demo script that demonstrates how to use the extended `interactionR()` function to perform competitive risk modeling with interaction terms.

## Requirements

- **R** version >= 4.0.0
- **cmprsk** package (for fitting competing risks regression models)

To install the required packages, run the following command in R:
```R
install.packages("cmprsk")
```

## Usage

1. **interactionR_crr.R**:
   - This script contains the core function that extends `interactionR()` to handle competitive risk models.
   - To use, source the script and run the function as per your dataset.

2. **interactionR_crr_demo.R**:
   - This script provides a step-by-step guide on how to apply the interactionR function on a sample dataset using the competitive risk regression (`crr`) model.
   - Simply source the script and run it to see an example of the function in action.

### Example Workflow

1. Source the main function:
```R
source("interactionR_crr.R")
```

2. Run the demo script to see how the function works:
```R
source("interactionR_crr_demo.R")
```

### Contribution

Feel free to contribute by submitting issues or pull requests. If you have suggestions or improvements, we encourage you to collaborate!

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
