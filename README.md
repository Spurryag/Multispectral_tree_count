# Application of Object Oriented Methods and Semi-Supervised Machine learning Algorithms for tree counts using Geospatial Data (MSc Thesis)

## Abstract

The recent democratisation of aerial image capture devices has led to an abundance of Geospatial data and created the need for innovative ways of analysing such data. The current analysis
methods are based around the use of specialised software to extract actionable intelligence from such datasets. While boasting a very high level of accuracy, the use of such software can
often require much time and computational power to provide insights for the end user. This thesis aimed to extract and classify different regions of a map based on how similar the spatial
points are between each other, in an automated fashion, through the use of the R Studio software suite. The study data encompassed the PortMoak Moss bog of Scotland, which covers an
area of 0.459 km squared. The objective of classifying the PortMoak Moss map regions relates to the identification of tree-covered and non-tree covered areas, such that an estimate of the
number of trees present is obtained. 

To do so, Object-oriented methods and a Semi-supervised machine learning framework were considered for the task. The results from the Edge detection methods, involving the use of a laplace edge detector and a sobel edge detector, indicate that a range of tree counts from 154 to 321 trees is present in the area. Conversely, the results from the semi-supervised
learning framework, including a K-means clustering algorithm and Random Forest classifier, suggest that around 3,684 trees are present. Nonetheless, given an initial estimated number of trees ranging from 436 to 9,687 based on certain assumptions and the area under study, the semi-supervised framework appears to provide a better estimation than the Object detection methods.
