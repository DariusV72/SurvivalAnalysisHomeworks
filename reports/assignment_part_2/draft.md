# Survival Analysis Project: Part II

## Report

### Simulation Algorithm

**Task**

Provide a step-by-step description of the simulation procedure.

For example:
(a) Generate covariates.
(b) Generate event times.
(c) Generate censoring times.
(d) Construct the observed dataset.
(e) Fit the proposed model.
(f) Store quantities of interest.
(g) Repeat for many simulation runs.

A flowchart is encouraged.

**Solution**

### Software Implementation

### Performance Measures

We evaluate the estimates of β_D, β_E, β_A and β_K, with β_E as the main target. In each of the R = 2000 runs we store the estimate β̂ and its standard error. With β̄ the average estimate over all runs, we compute:

Bias = β̄ − β: is the estimate systematically wrong?
Empirical SE = the standard deviation of the estimates: how much does the estimate vary between runs? (Variance = SE².)
MSE = the average of (β̂ − β)² ≈ Bias² + Variance: the total error.
Coverage = the share of 95% confidence intervals that contain the true β: are the standard errors reliable? It should be close to 95%.

For each scenario we also report the observed censoring rate and the number of runs that failed to converge (these are excluded).

Both the log-logistic and the Weibull fits are AFT models fitted with survreg, so their coefficients are on the same scale and can be compared with the same true values. When age is left out, only β_D, β_E and β_K are evaluated.

### Alternative Analysis

(a) Wrong distribution: fitting a Weibull model
We fit a Weibull AFT model instead of the log-logistic one.
Why it is wrong: the Weibull hazard can only increase, decrease or stay constant. It cannot rise and then fall, so it cannot capture the early "dangerous peak" in the true hazard.
(b) Omitting an important covariate: leaving out age
We fit the correct log-logistic model, but without age.
Why it is wrong: age is a confounder. It increases the number of comorbidities (K) and the chance of an extra-hepatic cause (E), and it also shortens survival directly.

### Hypotheses

> **Research questions:**
> RQ1. How accurately (bias, variance, MSE) does the log-logistic model estimate the true β ’s, and how does this change with sample size (n = 200 vs. 500) and loss to follow-up (5% vs. 10%)?
> RQ2. How biased are the estimates if we (a) fit a Weibull model, which cannot capture the early peak, or (b) leave out the confounder age?

#### RQ1: Correct log-logistic model

**H1 (bias).** When we fit the correct log-logistic model, the estimated effects will be very close to the true values in all four scenarios. Any small bias at n = 200 will become even smaller at n = 500. The reason is that the fitted model is the same as the model that generated the data, and the censoring does not depend on the survival time, so the estimates are correct on average.

**H2 (sample size).** The estimates will vary less from run to run with 500 patients than with 200. We expect the empirical standard error to drop by about one third, and the MSE to fall to roughly 40% of its value at n = 200. The reason is that more patients give more information. Since the bias is close to zero, the MSE is almost entirely caused by variance and falls with it.

**H3 (loss to follow-up).** Increasing loss to follow-up from 5% to 10% will make the estimates slightly less precise, but will not make them biased. This effect will be much smaller than the effect of sample size. The reason is that extra random censoring only removes some observed events. It does not push the estimates in a particular direction. Many patients are already censored at the 90-day cutoff, so losing an extra 5% makes little difference.

#### RQ2a: Wrong model (Weibull instead of log-logistic)

**H4.** The Weibull model will clearly misestimate the intercept and the shape of the survival curve. The covariate effects (days to treatment, age, comorbidities, etiology) will be only slightly biased. The reason is that the two models differ only in the assumed shape of the survival distribution. Most of the misfit is absorbed by the intercept and shape. The Weibull hazard can only go up or down over time, so it cannot reproduce the early peak in risk. Because many patients are censored, the covariate effects will still be somewhat affected.

**H5.** The bias of the Weibull model will not go away with more patients, while its variance will. This means that the Weibull model will look relatively worse at n = 500 than at n = 200, compared with the correct model. The reason is that a wrong model is precisely estimating the wrong thing: more data makes the estimate more stable but not more accurate.

**H6.** The Weibull model will be somewhat more biased at 10% loss to follow-up than at 5%. The reason is that the more of the data is censored, the more the estimates depend on the assumed shape of the distribution, which is wrong here.

#### RQ2b: Leaving out the confounder age

**H7.** Leaving age out of the model will make comorbidities and damage outside the liver look more harmful than they really are. The reason is that in our data, older patients have more comorbidities and more often have damage outside the liver, and older patients also survive shorter. Without age in the model, part of the harmful effect of age is wrongly attributed to comorbidities and etiology.

**H8.** The effect of days to treatment will stay close to its true value when age is left out, but it will be estimated slightly less precisely. The reason is that days to treatment was generated independently of age, so age does not confound it. The variation in survival caused by age becomes unexplained noise. This makes the estimates a bit less precise, but does not push them in a particular direction.


## Appendix

### Response to Feedback

> "The only comment I have is that you do not look up the literature after Gemini gave you information about TPE and liver failure. Can you check this for the next assignment?"

**Response:** As a fact-checking the provided information from Geimini Pro, based on ([https://www.aasld.org/liver-fellow-network/core-series/why-series/why-would-we-consider-plasma-exchange-acute-liver](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239959/#Sec12)), TPE is used for AFL for either as a "bridge" to liver recovery or maintain stability until liver transplant is possible. When it comes to comorbitidies, older age is strongly associated with an increased accumulation and severity of extrahepatic comorbid diseases (such as cardiovascular disease, chronic kidney disease, and type 2 diabetes) [Jepsen, 2014, pp. 7223–7225]. Metabolic comorbidities (obesity, type 2 diabetes, dyslipidemia) are the primary drivers of metabolic problems associated with steatotic liver disease [Huang et al., 2023, pp. 389–391]. Patients receiving TPE for acute liver failure are critically ill intensive-care patients and their trajectory is determined within days to weeks by either native hepatocyte recovery or emergency transplant, making outpatient hospital readmission an uninformative metric for acute TPE efficacy [Chris-Olaiya et al., 2021, pp. 905–906]. Reviewing this literature justified moving away from a hospital readmission framework and instead defining the failure event in Part I as transplant-free survival (time to death or liver transplant), observed from the first TPE cycle with right-censoring at 90 days [Larsen et al., 2016, pp. 70–72; Chris-Olaiya et al., 2021, p. 906].


### Use of Artificial Intelligence Tools

### Contribution Statement

### Rubric (not part of the submission)

- Feedback: 5
- Code: 1
- Correct implementation: 1
- Presentation of results: 1
- Interpretation and discussion: 1
- Conclusions and recommendations: 1
