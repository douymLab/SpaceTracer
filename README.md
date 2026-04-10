# SpaceTracer (NEW)
SpaceTracer is an open-source algorithm capable of accurately detecting mosaic SNVs, including both nuclear SNVs and mitochondria SNVs, directly from spatial transcriptomics data. 

![flowchart](./figures/flowchart.png)

## Release Notes
- 2025/03/31: Version 1.0.0  
This is the initial version of SpaceTracer.
- 2026/02/25: Version 1.1.0  
This release focuses on updating the genotype calculation and enhancing the features used in the random forest model. Additionally, we've added more filter steps to improve the accuracy of the results.

We expect to release the **next version** of SpaceTracer in roughly one month. This release will feature substantial performance improvements.

**Key Improvements for Upcoming Release**
- Integrated Lysis Error Calculation \
Lysis error calculation for single samples will be incorporated directly into the full algorithm pipeline, providing more comprehensive and streamlined analysis.
- One-Command Execution Mode \
A new simplified execution option will allow users to run SpaceTracer with a single command, removing the dependency on Snakemake and making the workflow more accessible.
- Up to 10× Speed Improvement \
Major performance optimizations will enable processing speeds up to ten times faster than the previous version.


## Contact:
If you have any questions please contact us:  
Zhirui Yang: yangzhirui@westlake.edu.cn  
Mengdie Yao: yaomengdie@westlake.edu.cn  
Yanmei Dou: douyanmei@westlake.edu.cn
