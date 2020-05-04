# Predicting effective reproduction number, Rt, for Pakistan
**By Faraz Ahmed Khan and Musabbir Abdul Majeed**

We implemented the renewable model framework proposed by Cori et al.[1] and extended by Parag et al.[2]. This approach allows us to estimate instantaneous reproduction number from the past incidence or infection curve, and predict infections for the next day.

Our preliminary results for Balochistan, KPK, Punjab, and Sindh are shown below and they will be updated daily. In the figures, the top and bottom graphs on each page show the instantaneous reproduction number and the number of new infections respectively. The blue dots are the number of new infections reported by the Pakistan government, http://covid.gov.pk/stats/pakistan. The black curve is the mean estimate and the grey shaded region is the 95% confidence interval.

To verify the applicability of the renewable model framework in the context of Pakistan and our hyperparameters, we have predicted the number of infections for the last six days (grey region) and compared with the reported number of infections.

[1. A new framework and software to estimate time-varying reproduction numbers during epidemics](https://academic.oup.com/aje/article/178/9/1505/89262)

[2. Optimising Renewal Models for Real-Time Epidemic Prediction and Estimation](https://www.biorxiv.org/content/10.1101/835181v1.abstract)

### Effective reproduction number of Balochistan

![Effective reproduction number of Balochistan](/pakistan_data/Predictions/Balochistan.png)

### Effective reproduction number of KPK

![Effective reproduction number of KPK](/pakistan_data/Predictions/KPK.png)

### Effective reproduction number of Punjab

![Effective reproduction number of Punjab](/pakistan_data/Predictions/Punjab.png)

### Effective reproduction number of Sindh

![Effective reproduction number of Sindh](/pakistan_data/Predictions/Sindh.png)