# Compare-DiD-Estimators

**Case**: staggered treatments, homogenous (constant) TEs. (100 units, 15 time periods)
- "Cohort 1": 35 units are treated once at t=7,  constant ATT(g=1,K)=49 for K=0-8 periods after the treatment (t=7-15)
- "Cohort 2": 35 units are treated once at t=11, constant ATT(g=2,K)=49 for K=0-4 periods after the treatment (t=11-15)
- 30 units are never treated
![image](https://user-images.githubusercontent.com/7883904/154598476-5d299307-df4b-4d43-a16e-8d2e41afa93e.png)
![image](https://user-images.githubusercontent.com/7883904/154598584-d81ba073-ae4d-40f0-b43e-5bd4930679b2.png)


**Case**: staggered treatments, heterogenous (constant) TEs. (100 units, 15 time periods)
- "Cohort 1": 35 units are treated once at t=7,  constant ATT(g=1,K)=20 for K=0-8 periods after the treatment (t=7-15)
- "Cohort 2": 35 units are treated once at t=11, constant ATT(g=2,K)=40 for K=0-4 periods after the treatment (t=11-15)
- 30 units are never treated
![image](https://user-images.githubusercontent.com/7883904/154596389-2156dde9-ea30-4a67-b57e-4b11ec6fe271.png)
![image](https://user-images.githubusercontent.com/7883904/154597194-c0cdf7e7-8497-4068-8655-cd8e14e8068b.png)
- TE estimates for K>4 periods post-treatment are non-sense b/c for the 2nd cohort, we only have data up to 4-periods post treatment, so estimates of ATT_K for K>4 only average out the TE for the first cohort (which has a smaller ATT=20)
- TWFE-OLS & CDLZ are biased. 
