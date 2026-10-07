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



### Hypotheses

> **Research questions:**
> RQ1. How accurately (bias, variance, MSE) does the log-logistic model estimate the true β ’s, and how does this change with sample size (n = 200 vs. 500) and loss to follow-up (5% vs. 10%)?
> RQ2. How biased are the estimates if we (a) fit a Weibull model, which cannot capture the early peak, or (b) leave out the confounder age?

## Appendix

### Response to Feedback

> "The only comment I have is that you do not look up the literature after Gemini gave you information about TPE and liver failure. Can you check this for the next assignment?"

**Response:** As a fact-checking the provided information from Geimini Pro, based on ([https://www.aasld.org/liver-fellow-network/core-series/why-series/why-would-we-consider-plasma-exchange-acute-liver](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239959/#Sec12)), TPE is used for AFL for either as a "bridge" to liver recovery or maintain stability until liver transplant is possible. When it comes to comorbitidies, older age is strongly associated with an increased accumulation and severity of extrahepatic comorbid diseases (such as cardiovascular disease, chronic kidney disease, and type 2 diabetes) [Jepsen, 2014, pp. 7223–7225]. Metabolic comorbidities (obesity, type 2 diabetes, dyslipidemia) are the primary drivers of metabolic problems associated with steatotic liver disease [Huang et al., 2023, pp. 389–391]. Patients receiving TPE for acute liver failure are critically ill intensive-care patients and their trajectory is determined within days to weeks by either native hepatocyte recovery or emergency transplant, making outpatient hospital readmission an uninformative metric for acute TPE efficacy [Chris-Olaiya et al., 2021, pp. 905–906]. Reviewing this literature justified moving away from a hospital readmission framework and instead defining the failure event in Part I as transplant-free survival (time to death or liver transplant), observed from the first TPE cycle with right-censoring at 90 days [Larsen et al., 2016, pp. 70–72; Chris-Olaiya et al., 2021, p. 906].


### Contribution Statement

### Rubric (not part of the submission)

- Feedback: 5
- Code: 1
- Correct implementation: 1
- Presentation of results: 1
- Interpretation and discussion: 1
- Conclusions and recommendations: 1

**Initial answer for feedback:**

**Response:** In Part I, we utilized exploratory AI prompts to brainstorm clinical contexts without subsequently verifying the outputs against published medical literature. Below, we fact-check each of the three exploratory questions against peer-reviewed hepatology literature, grounding our simulation design, covariate selections, and survival event definitions:

#### Question 1: In what cases is TPE used for liver failure? Is this a permanent treatment/cure or does it just postpone the inevitable?

##### Key Reference

- **Chris-Olaiya, A., Kapoor, A., Ricci, K. S., & Lindenmeyer, C. C. (2021).** Therapeutic plasma exchange in liver failure. _World Journal of Hepatology_, 13(8), 904–915. [DOI: 10.4254/wjh.v13.i8.904](https://doi.org/10.4254/wjh.v13.i8.904) (PMID: 34552697; PMCID: PMC8422921)

##### How it answers the question

- **Established Clinical Indications:**
  - Acute Liver Failure (ALF): High-volume TPE holds an ASFA Category I (Grade 1A) indication across drug-induced (e.g., acetaminophen toxicity), viral, toxic, and autoimmune etiologies [Chris-Olaiya et al., 2021, Table 1 & pp. 905–907].
  - Wilson's Disease: TPE rapidly clears toxic circulating free copper, abating severe intravascular hemolysis and secondary renal failure [Chris-Olaiya et al., 2021, p. 907].
  - Acute-on-Chronic Liver Failure (ACLF): TPE is applied as extracorporeal liver support (ASFA Category III) to clear inflammatory cytokines, bilirubin, and ammonia, and to improve coagulopathy [Chris-Olaiya et al., 2021, pp. 908–910].
- **Cure vs. Postponing the Inevitable (Bridge):**
  - Not a permanent structural cure: TPE removes toxins and inflammatory mediators from circulation, but it cannot reverse established underlying cirrhosis or permanent parenchymal damage [Chris-Olaiya et al., 2021, p. 905].
  - ALF as a Bridge to Recovery: In acute liver failure with previously healthy hepatic tissue, dampening systemic inflammation allows native hepatocytes to regenerate, achieving true transplant-free survival (58.7% with TPE vs. 47.8% with standard medical therapy in Larsen et al., 2016) [Chris-Olaiya et al., 2021, p. 906].
  - ACLF / Cirrhosis as a Bridge to Transplant: In end-stage chronic liver failure, TPE acts strictly as a temporary bridge to emergency liver transplantation by stabilizing multi-organ failure; without a transplant or resolution of the precipitating insult, long-term mortality remains high [Chris-Olaiya et al., 2021, pp. 908–911].

#### Question 2: How does age affect comorbidities, how do comorbidities affect etiology, and vice versa?

##### Key Reference

- **Jepsen, P. (2014).** Comorbidity in cirrhosis. _World Journal of Gastroenterology_, 20(23), 7223–7230. [DOI: 10.3748/wjg.v20.i23.7223](https://doi.org/10.3748/wjg.v20.i23.7223) (PMID: 24966593; PMCID: PMC4064072)
- _(Supplementary on etiology trends)_ **Huang, D. Q., Terrault, N. A., Tacke, F., et al. (2023).** Global epidemiology of cirrhosis — aetiology, trends and predictions. _Nature Reviews Gastroenterology & Hepatology_, 20(6), 388–398. [DOI: 10.1038/s41575-023-00759-2](https://doi.org/10.1038/s41575-023-00759-2)

##### How it answers the question

- **Age -> Comorbidities:**
  - Older age is strongly associated with an increased accumulation and severity of extrahepatic comorbid diseases (such as cardiovascular disease, chronic kidney disease, and type 2 diabetes) [Jepsen, 2014, pp. 7223–7225].
  - In survival models, age is a major confounder because it directly drives comorbidity burden while independently increasing overall mortality risk (${K \leftarrow A \rightarrow T}$) [Jepsen, 2014, pp. 7224–7226].
- **Comorbidities -> Etiology & Liver Damage:**
  - Metabolic comorbidities (obesity, type 2 diabetes, dyslipidemia) are primary causative drivers of metabolic dysfunction-associated steatotic liver disease (MASLD/MASH) [Huang et al., 2023, pp. 389–391].
  - Concurrent metabolic and cardiovascular comorbidities accelerate liver fibrosis and increase mortality across both viral and alcohol-associated liver disease etiologies [Jepsen, 2014, p. 7225; Huang et al., 2023, p. 392].
- **Etiology -> Comorbidities:**
  - Underlying liver disease etiology directly shapes patient comorbidity patterns: alcohol-associated liver disease is frequently coupled with malnutrition, cardiomyopathy, and chronic pancreatitis, whereas chronic viral hepatitis is tied to specific extrahepatic autoimmune manifestations [Jepsen, 2014, pp. 7225–7227].

#### Question 3: Using hospital readmissions as a study example: What is the role and impact of plasma exchange in patients with liver failure?

##### Key Reference

- **Chris-Olaiya, A., Kapoor, A., Ricci, K. S., & Lindenmeyer, C. C. (2021).** Therapeutic plasma exchange in liver failure. _World Journal of Hepatology_, 13(8), 904–915. [DOI: 10.4254/wjh.v13.i8.904](https://doi.org/10.4254/wjh.v13.i8.904) (PMID: 34552697; PMCID: PMC8422921)
- _(Landmark trial)_ **Larsen, F. S., et al. (2016).** High-volume plasma exchange in patients with acute liver failure: An open randomised controlled trial. _Journal of Hepatology_, 64(1), 69–78. [DOI: 10.1016/j.jhep.2015.08.018](https://doi.org/10.1016/j.jhep.2015.08.018)

##### How it answers the question

- **Primary Clinical Endpoints for TPE in Liver Failure:**
  - Clinical trials and guidelines for TPE in acute liver failure focus on **transplant-free survival, 30/90-day survival, and bridging to emergency liver transplantation**, rather than post-discharge hospital readmission [Larsen et al., 2016, pp. 69–71; Chris-Olaiya et al., 2021, pp. 905–907].
- **Why Readmission is Not a Primary Outcome for TPE:**
  - Patients receiving TPE for acute liver failure are critically ill intensive-care patients; their clinical trajectory is determined acutely within days to weeks by either native hepatocyte recovery or fatal multi-organ failure / emergency transplantation, making outpatient hospital readmission an uninformative metric for acute TPE efficacy [Chris-Olaiya et al., 2021, pp. 905–906].
- **Justification for Our Survival Analysis Design:**
  - Reviewing this literature justified moving away from a hospital readmission framework and instead defining the failure event in Part I as **transplant-free survival (time to death or liver transplantation)**, observed from the first TPE cycle with right-censoring at 90 days [Larsen et al., 2016, pp. 70–72; Chris-Olaiya et al., 2021, p. 906].

### Use of Artificial Intelligence Tools
