# hierarchical_clustering
My own implementation of the bottom-up hierarchical clustering algorithm.

This program runs the algorithm until a specified number, k, of clusters remain and output a description of each of these k clusters. 

It takes as input a small gene expression dataset.

The algorithm uses Euclidean distance as the distance between two expression profiles, and is able to do single-link, average-link, or complete-link when determining the distances between clusters.

The expression data file contains a table in Tab Separated Values (tsv) format. Each line of the file describes a single gene, and individual columns inside each line, separated by tabs, depict the different aspects of the gene.

The first column lists the identifier of the gene. The second column lists a common name for the gene and a description of its function. The remaining columns list measurements for the gene under various conditions. These measurements are log expression ratios. Note that this program only uses these measurement values in doing the clustering.

A description of the k clusters is printed by the program. The description of a cluster lists, by lines and separated by tabs, the identifier, the name and description, and the average value of the measurements for genes contained in that cluster. It also outputs in a new line the average value of all the measurements in that cluster. 

To ensure the uniqueness of the output, the program orders clusters by the average measurement values, in ascending order, and within each cluster, genes are also ordered by their average measurement value, in ascending order.

To run the program:

Rscript cluster.R test_cases/tiny-yeast.tsv link_type k

Where link_type is one of the following: single-link (S), average-link (A), or complete-link (C). k is the number of clusters to return.

