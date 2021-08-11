RATS Library
===========

**R**andomness **A**uto-**T**est **S**uite (**RATS**) is a proof-of-concept implementation of the statistical test for random and pseudo-random number generators for cryptographic applications.

The following schemes are currently supported...

1. Frequency (Monobit) Test
2. Frequency Test within a Block 
3. Runs Test
4. Tests for the Longest-Run-of-Ones in a Block
5. Binary Matrix Rank Test
6. Discrete Fourier Transform (Spectral) Test
7. Non-overlapping Template Matching Test
8. Overlapping Template Matching Test
9. Maurer's "Universal Statistical" Test
10. Linear Complexity Test
11. Serial Test
12. Approximate Entropy Test
13. Cumulative Sums (Cusums) Test
14. Random Excursions Test
15. Random Excursions Variant Test 



# Warning

This code is not to be considered secure, efficient or fully portable. 
Its purpose is not to be used for practical engineering applications, but to provide a tool for the research community to verify, analyze and reproduce the implementation efficiency.

# How to use?

# Copyright

Copyright (C) Wuhan University 2019


# Source Control
"all" 更新了cuda实现的LinearComplexity，采用并行优化，具体优化在于分块，并行处理各个块的线性复杂度。

    sts-2.1。2中，进行1000轮1000000bit测试使用的时间为5005.1961s。
    “all” 更新之后，使用的时间为1387.165s。

"add nonOverlappingTemplateMatchings" 

    在文件nonOverlappingTemplateMatchings.c中添加 N O N O V E R L A P P I N G  T E M P L A T E  T E S T ，还未做优化，该测试在进行1000000bit数据量的情况下，耗时约2.103s。

"2021.8.11 1:20"

    添加了OverlappingTemplateMatchings.c randomExcursions.c randomExcursionsVariant.c三种测试方法
    在未做优化的情况下，做完全的19组实验，使用1000000bit数据，所花费的时间为3.6145s，对比sts：5.3788s

"2021.8.11 14:31"

    为nonOverlappingTemplateMatchings进行了并行优化，按照分块并行计算，其余地方并没有进行修改，只是将最终p-value的计算结果移至cpu进行操作，
    在gpu中只能进行简单的运算。
