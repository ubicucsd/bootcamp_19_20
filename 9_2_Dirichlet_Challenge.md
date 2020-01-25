## Dirichlet Process Means Clustering

Copy the starter file from 
```/home/ubuntu/clustering_lesson/dmp.py``` 

The dirichlet process is just a method of clustering without knowing the number of clusters ahead of time. Something like this algorithm may be used for example 2 in the list of clustering related bioinformatics problems. Here are the steps:
 
&nbsp;&nbsp;  I. Add the first point in your data as the first center
  
&nbsp;&nbsp;  II. Iterate through all of the points in your data, calculating the distance between the current point and all existing centers. 
  
&nbsp;&nbsp;&nbsp;&nbsp; A. If the point falls within a certain threshold distance of the center, add that point to that center's cluster
  
&nbsp;&nbsp;&nbsp;&nbsp; B. If the point does not fall within a certain threshold distance of the center, add that point to the list of centers
    
&nbsp;&nbsp; III. Once you have gone through the Dirichlet Process, you do the means part. Redifine the centers of each cluster to be the average of all points in the cluster, like you would in most standard clustering algorithms. 
  
&nbsp;&nbsp; IV. Repeat steps II and III either until an iteration provides no change in the assignment of points to centers, or a maximum number of cycles is reached. 

**HINT:** I recommend using arrays to store the cluster information! You will need an array that is the same length as the number of points(index 1 is point 1 and so on). Then, you will assign a cluster number to each index as you go along (keep track of the number of clusters//increment when you create a new one). Also, you need an array of the actual centers. 
  
