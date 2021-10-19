## Two dimensional hydrodynamic model

### Explanations of model 

For the two dimensional hydrodynamic models, the forces on netting are calculated based on individual a virtual net panel of netting. The twines and knots in the net panel are considered as an integrated structure. In this module, the net panel is defined by three nodes because three (non-collinear) points can determine a unique plane in Euclidean geometry. In practice, the force is usually decomposed into two components: drag force ![formula](https://render.githubusercontent.com/render/math?math=\vec{F_D}) and lift force ![formula](https://render.githubusercontent.com/render/math?math=\vec{F_L}) (Cheng *et al.*, 2020).

![Fig.4](./figures/Fig.4.png)


![formula](https://render.githubusercontent.com/render/math?math=\quad\Large\vec{F_D}=0.5C_{D}\rho_{w}A_{t}\left|\vec{U}-\vec{v}\right|^{2}\vec{i_D})


![formula](https://render.githubusercontent.com/render/math?math=\quad\Large\vec{F_L}=0.5C_{L}\rho_{w}A_{t}\left|\vec{U}-\vec{v}\right|^{2}\vec{i_L})


* ![formula](https://render.githubusercontent.com/render/math?math=\rho_{w}) is the fluid density.

* ![formula](https://render.githubusercontent.com/render/math?math=A_t) is the area of a virtual net panel (*i.e*., the area of the triangular P1-P2-P3 in the above figure).

* ![formula](https://render.githubusercontent.com/render/math?math=\vec{U}) is the undisturbed incoming flow velocity in the upstream of the net panel.

* ![formula](https://render.githubusercontent.com/render/math?math=\vec{v}) is the velocity of the structure.

* The unit vectors ![formula](https://render.githubusercontent.com/render/math?math=\vec{i_D}) and ![formula](https://render.githubusercontent.com/render/math?math=\vec{i_L}) which are used to indicate the directions of forces are defined by the following equations (![formula](https://render.githubusercontent.com/render/math?math=\vec{e_n}) is the unit normal vector of the virtual net panel):


![formula](https://render.githubusercontent.com/render/math?math=\Large\vec{i_D}=\frac{\vec{U}-\vec{v}}{\left|\vec{U}-\vec{v}\right|})

![formula](https://render.githubusercontent.com/render/math?math=\Large\vec{i_L}=\frac{(\vec{U}-\vec{v})\times\vec{e_n}\times(\vec{U}-\vec{v})}{\left|(\vec{U}-\vec{v})\times\vec{e_n}\times(\vec{U}-\vec{v})\right|})




* ![formula](https://render.githubusercontent.com/render/math?math=\vec{C_D}) and ![formula](https://render.githubusercontent.com/render/math?math=\vec{C_L}) are the drag and lift force coefficients in the two dimensional models, respectively. These force coefficients are usually obtained from experiments that approximate the ideal conditions of a finite net panel in an infinite flow field.

### How to implement in code

Example code:

```python
import hydroModules as hdm

netting = hdm.two_dimensional.netting("S3",
                                      meshinfo['surfs_netting'],
                                      caseinfo['Net']['Sn'],
                                      caseinfo['Net']['twineDiameter'],
                                      )
                                      
```

> **model_index**: [string] [-]. To indicate the model function, *e.g*., 'S1', 'S2', 'S3'.
>
> **hydro_element**: [list] [-]. A python list to indicate how the net panel are connected. *e.g.*, [[p1,p2,p3][p2,p3,p4,p5]...]. If the input net panel contains 4 nodes, it will automaticly decomposed to 3 node net panel.
>
> **solidity**: [float] [-]. The solidity of netting.
>
> **dw0**: [float] [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
