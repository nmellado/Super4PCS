# Datasets {#datasets}

If you want to go further, you can also download the demo datasets, including 3D models from the original paper. To do so, you can either call `cmake dl-datasets` if you use Super4PCS from sources, or directly at downloading this [file](https://www.irit.fr/~Nicolas.Mellado/dl/datasets/point-clouds/Super4PCS-dataset-demo1.zip) (MD5: `ad1e172902b41a3f17e9b4901adf3ba5`). To run the demo,
 * Go to `assets/demo1/`,
 * Call the `run` script, which will register the 3D models (point clouds, meshes, etc).
 * To check the registration output, you can use the meshlab projects (files `*.mlp`) provided in the dataset subfolders.
   * Note that this demo does not include local registration, ICP must be ran in post-process to obtain fine alignment.
