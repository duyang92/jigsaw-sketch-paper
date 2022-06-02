# Jigsaw-Sketch

Jigsaw-Sketch is a new algorithm for finding the top-k elephant flows. It is accepted by Science China Information Sciences. This is the repository for the source codes.

## Introduction
Finding top-$k$ elephant flows in high-speed networks is one of the most fundamental network measurement tasks. It is more challenging than per-flow size estimation since the IDs and sizes of top-$k$ flows must be tracked simultaneously. Most existing studies only record the IDs of a small number of elephant flows to fit their estimators in the extremely limited high-speed on-chip memory. However, these solutions need too many memory accesses when a packet arrives to track the elephant flows with high accuracy, which limits their practicability. Therefore, this paper proposes Jigsaw-Sketch, a new algorithm to find the top-$k$ elephant flows with much fewer memory accesses while achieving high memory efficiency and accuracy. In this design, we propose a novel two-stage jigsaw storage scheme, which can capture the candidate top-$k$ flows from massive network steams efficiently, and further find the top-$k$ elephant flows with high memory efficiency and only a few memory accesses for each packet. Extensive experimental results based on real network traces show that Jigsaw-Sketch improves the packet processing throughput by at least $86\%$, while achieving smaller memory footprints and higher accuracy compared to the SOTA.

## About this repo
On the CPU platform, the codes of Jigsaw-Sketch and the baselines are compiled using gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04). We employ the default optimization options when compiling codes. For each algorithm, the compiling command is as follows.
```shell
g++ main.cpp MurmurHash3.cpp -O -m64
```

