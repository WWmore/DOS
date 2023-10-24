# Discrete Orthogonal Structures 

Felix Dellinger, Xinye Li, and Hui Wang* (corresponding author)<br>
| [Webpage](https://www.huiwang.me/projects/10_project/) | [Full Paper](https://www.huiwang.me/assets/pdf/2023SMI.pdf) | [Publication Page](https://doi.org/10.1016/j.cag.2023.05.024) |<br>
![Teaser image](assets/teaser.png)

Abstract: *To represent smooth geometric shapes by coarse polygonal meshes, visible edges often follow special families of curves on a surface to achieve visually pleasing results. Important examples of such families are principal curvature lines, asymptotic lines or geodesics. In a surprisingly big amount of use-cases, these curves form an orthogonal net. While the condition of orthogonality between smooth curves on a surface is straightforward, the discrete counterpart, namely orthogonal quad meshes, is not. In this paper, we study the definition of discrete orthogonality based on equal diagonal lengths in every quadrilateral. We embed this definition in the theory of discrete differential geometry and highlight its benefits for practical applications. We demonstrate the versatility of this approach by combining discrete orthogonality with other classical constraints known from discrete differential geometry. Orthogonal multi-nets, i.e. meshes where discrete orthogonality holds on any parameter rectangle, receive an in-depth analysis.*

This repository contains the implementation associated with the paper "Discrete Orthogonal Structures", which can be found [here](https://doi.org/10.1016/j.cag.2023.05.024). 
Please cite the paper if you use this code in your project. 

<section class="section" id="BibTeX">
  <div class="container is-max-desktop content">
    <h2 class="title">BibTeX</h2>
    <pre><code>@Article{DOS2023,
      author       = {Dellinger, Felix and Li, Xinye and Wang, Hui},
      title        = {Discrete Orthogonal Structures},
      journal      = {Computers & Graphics},
      volume       = {114},
      pages        = {126--137},
      month        = {June},
      year         = {2023},
      doi          = {10.1016/j.cag.2023.05.024},
      url          = {https://www.huiwang.me/projects/10_project/}
}</code></pre>
  </div>
</section>


## Set up a working environment in Windows / MacOS

Using Anaconda to install every package.

    1. Download Anaconda

    2. Open Anaconda Prompt
    ```
    $ conda create -n geo 
    $ conda activate geo
    $ conda install mayavi traits traitsui qt pyqt vtk scipy spyder 
    $ conda install -c haasad pypardiso
    ```
    3. Open Anaconda, under "geo" environment open Spyder

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


## File Relations

![File](assets/files.png)

<details>
<summary><span style="font-weight: bold;">File tree.</span></summary>

  ![File](assets/tree.png)

</details>
<br>

![File](assets/frame.png)

<details>
<summary><span style="font-weight: bold;">File notice.</span></summary>

  - note: in the above figure, rectangles represent .py file and rounded rectangles represent folder names.

  - files in geometrylab folder are basic, nothing need to be changed.

  - if you want to test how it works, just try python files in geometrylab/test: ex. run paneling.py, then a GUI window will be opened.

  - if you want to add a new optimization project, you can create a new folder similar to the folder 'archgeolab'. Then the mesh geometry, optimization and GUI will be based on the files in geometrylab folder (please refer to archgeolab).

  - archgeolab/archgeometry: meshpy.py --> quadrings.py --> gridshell_new.py --> gui_basic.py --> guidedprojection_orthonet.py --> opt_gui_orthonet.py --> readfile_orthonet.py

</details>
<br>


## Mesh Optimization
The optimizer uses Guided Projection Algorithm, a Gauss-Newton algorithm, as dissused in the paper [Form-finding with polyhedral meshes made simple](https://doi.org/10.1145/2601097.2601213), in a Python environment to produce quadmesh models.

<details>
<summary><span style="font-weight: bold;">Abstract of the paper 'Form-finding with Polyhedral Meshes Made Simple'</span></summary>

  We solve the form-finding problem for polyhedral meshes in a way which combines form, function and fabrication; taking care of user-specified constraints like boundary interpolation, planarity of faces, statics, panel size and shape, enclosed volume, and last, but not least, cost. Our main application is the interactive modeling of meshes for architectural and industrial design. Our approach can be described as guided exploration of the constraint space whose algebraic structure is simplified by introducing auxiliary variables and ensuring that constraints are at most quadratic. Computationally, we perform a projection onto the constraint space which is biased towards low values of an energy which expresses desirable "soft" properties like fairness. We have created a tool which elegantly handles difficult tasks, such as taking boundary-alignment of polyhedral meshes into account, planarization, fairing under planarity side conditions, handling hybrid meshes, and extending the treatment of static equilibrium to shapes which possess overhanging parts.

</details>
<br>

## ArchGeo Visualization


Python programming visualization for optimization problems in the Architectural Geometry / Geometry Processing area.
This implementation major works on quad meshes.


![File](assets/mayavi.png)


### Implementation layout
[![layout](assets/layout.png)](https://www.youtube.com/embed/1l6DCW9BmYM)


### Implementation of a principal net optimized from an orthogonal PQ mesh
[![PQ](assets/pq.png)](https://www.youtube.com/embed/m-CFC0XZ488)


### Implementation of a minimal net optimized from an orthogonal A-net
[![Anet](assets/anet.png)](https://www.youtube.com/embed/KQbJ2e_Ow7M)


### Implementation of a CMC net optimized from an orthogonal S-net with const. radius
[![CMC](assets/cmc.png)](https://www.youtube.com/embed/vgb9A6uAidw)


### Implementation of a principal stress net from an orthogonal equilibrium mesh
[![Funicular](assets/funicular.png)](https://www.youtube.com/embed/sOzjRHIrR-s)


------
## Contributions
If you find this codebase and paper helpful in your research, welcome to cite the paper and give a :star: .
This project was initially developed by [Davide Pellis (davidepellis@gmail.com)](https://scholar.google.com/citations?user=JnocFM4AAAAJ&hl=en). 
It is still under development. 
Please feel free to push issues or submit requests to contribute to our codebase (e.g. add shadows of the model).
Welcome to work together to make a Grasshopper plugin if you are interested.
Hoping this codebase is helpful for your research work. 

Welcome to the research collaborations.
