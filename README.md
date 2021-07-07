# nls_for_EPONs
This project shows how to build models for delay percentiles in EPONs, using machine learning techniques 
for the accurate characterization of delay percentiles.

There are 6 datasets (csv files) for Poisson traffic where we have considered 1- and 10-Gbps EPON comprising eight ONUs. 
We have examined EPONs with different fiber lengths being set to 4, 8, and 20 km. 
The guard band is set to 128 ns. 
To arbitrate the channel access in the upstream, we have implemented the well-known interleaved 
polling with adaptive cycle time (IPACT) protocol. 

The incoming traffic follows a Poisson arrival model with various traffic intensities 
(load varying between 0.2 and 0.99), while packet sizes follow the classical trimodal 
packet distribution observed at the  Amsterdam Internet Exchange (AMS-IX); namely 
packet sizes of 40, 576, and 1500 Bytes with percentages of 7/12, 4/12, and 1/12, respectively. 

In total, the experiments include six datasets, three at 1 Gb/s (4km, 8km and 20 km), and 
another three at 10 Gb/s (again 4km, 8km and 20 km). The former three datasets 
include 10,000 packet delay measurements per dataset, while the later three contain 50,000 packet delay 
measurements per dataset. We note that these numbers were deemed sufficient to reach the desired statistical stability.

With such individual packet delay values, we construct the ML models both for the average delay and also 
for different delay percentiles, as explained in the next sections. Delay percentiles are useful when the 
goal is to dimension PON scenarios where certain worse-case delay guarantees are provided. For instance, 
the 90th delay percentile represents the delay value experienced by the top-highest 10% delay.

There is also a chaoticmodel file to generate LRD traffic following the chaotic map model, and multiple LRD traces of PON delays generated under LRD traffic.

