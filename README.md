# Discrete Orthogonal Structures 

Felix Dellinger, Xinye Li, and Hui Wang* (corresponding author)<br>
| [Webpage](https://www.huiwang.me/projects/10_project/) | [Full Paper](https://www.huiwang.me/assets/pdf/2023SMI.pdf) | [Publication Page](https://doi.org/10.1016/j.cag.2023.05.024) |<br>
![Teaser image](assets/teaser.png)

This repository contains the implementation associated with the paper "Discrete Orthogonal Structures", which can be found [here](https://doi.org/10.1016/j.cag.2023.05.024). 

Abstract: *To represent smooth geometric shapes by coarse polygonal meshes, visible edges often follow special families of curves on a surface to achieve visually pleasing results. Important examples of such families are principal curvature lines, asymptotic lines or geodesics. In a surprisingly big amount of use-cases, these curves form an orthogonal net. While the condition of orthogonality between smooth curves on a surface is straightforward, the discrete counterpart, namely orthogonal quad meshes, is not. In this paper, we study the definition of discrete orthogonality based on equal diagonal lengths in every quadrilateral. We embed this definition in the theory of discrete differential geometry and highlight its benefits for practical applications. We demonstrate the versatility of this approach by combining discrete orthogonality with other classical constraints known from discrete differential geometry. Orthogonal multi-nets, i.e. meshes where discrete orthogonality holds on any parameter rectangle, receive an in-depth analysis.*

<section class="section" id="BibTeX">
  <div class="container is-max-desktop content">
    <h2 class="title">BibTeX</h2>
    <pre><code>@Article{SMI2023,
      author       = {Dellinger, Felix and Li, Xinye and Wang, Hui},
      title        = {Discrete Orthogonal Structures},
      journal      = {Computers & Graphics},
      volume       = {114},
      pages        = {126--137},
      month        = {June},
      year         = {2023},
      url          = {https://www.huiwang.me/projects/10_project/}
}</code></pre>
  </div>
</section>




## Instructions to set up a working environment in Windows / MacOS

Using Anaconda to install every package.

    1. Download Anaconda

    2. Open Anaconda Prompt in Windows searching
    ```
    $ conda create -n geo 
    $ conda activate geo
    $ conda install mayavi traits traitsui qt pyqt vtk scipy spyder 
    $ conda install -c haasad pypardiso
    ```
    3. Open Anaconda, under geo environment open Spyder

Once above installation failed because of versions conflict, then try below installations:

    ```
    $ conda create -n geo python=3.6
    $ conda activate geo
    $ pip install numpy scipy
    $ pip install python-vtk
    $ pip install mayavi --no-cache
    $ conda install -c haasad pypardiso
    $ conda install pyface
    ```

Or use the exported files within ```./conda/``` to set your environment

    ```
    $ conda env create -f environment.yml
    ```


<details>
<summary><span style="font-weight: bold;">Files instruction.</span></summary>
    1. files in geometrylab folder are basic, nothing need to be changed.

    2. if want to test how it works, try python files in geometrylab/test: ex. run paneling.py, then a GUI window will be opened.

    3. if want to add a new project, create a new folder named 'archgeolab'. the mesh geometry, optimization and GUI will be based on the files in geometrylab.

    4. archgeolab/archgeometry: meshpy.py --> quadrings.py --> gridshell.py --> gui_basic.py --> project folder (proj_orthonet)

    5. archgeolab/proj_orthonet: guidedprojection_orthonet.py --> opt_gui_orthonet.py --> readfile_orthonet.py
</details>
<br>