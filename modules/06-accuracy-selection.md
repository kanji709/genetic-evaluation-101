# Module 06 — Accuracy, selection and reporting

## Learning objectives

- Compute EBV accuracy and reliability from prediction error variance
- Predict response to selection with the breeder's equation
- Explain the trade-off between genetic gain and inbreeding
- Communicate EBVs clearly to breeders and non-specialists

## Core concepts

**Accuracy and reliability.**

  r = √(1 − PEV / σ²a),  reliability = r²

PEV comes from the inverse of the MME coefficient matrix (see Module 03 script, Step 4). For inbred animals, replace σ²a with (1 + F)σ²a. More records, more relatives with records and higher h² all raise accuracy.

**Breeder's equation.**

  R = i × r × σA / L

- i: selection intensity (how strictly you select)
- r: accuracy of selection
- σA: additive genetic standard deviation
- L: generation interval (years)

Genomics increases genetic gain mainly by raising r for young animals and shortening L.

**Inbreeding and diversity.** Selecting only the top EBVs concentrates on few families. Tools such as optimal contribution selection balance gain against the rate of inbreeding (ΔF).

**Reporting EBVs.** EBVs are relative to a genetic base, expressed in trait units, and should always come with accuracy. Breeders often see them combined into an economic **selection index**.

## Worked example

From the Module 03 output, rank the animals by EBV and compare accuracies. Then estimate annual response for a trait with σA = 9.5 kg, i = 1.4 (top ~20% selected), r = 0.6 and L = 3 years:

  R = 1.4 × 0.6 × 9.5 / 3 ≈ 2.7 kg per year

## Exercises

1. If genomic selection raises r to 0.75 and lowers L to 2 years, what is the new annual response?
2. Two bulls have EBVs of +12 (r = 0.45) and +10 (r = 0.95). Which would you recommend to a risk-averse farmer, and why?
3. Draft a one-paragraph explanation of "EBV" for a farmer newsletter.

## Further reading

- Rendel & Robertson (1950) Estimation of genetic gain in milk yield by selection in a closed herd of dairy cattle. *J. Genetics* 50: 1–8
- Meuwissen (1997) Maximizing the response of selection with a predefined rate of inbreeding. *J. Anim. Sci.* 75: 934–940
