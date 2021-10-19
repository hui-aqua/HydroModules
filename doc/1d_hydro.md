## One dimensional hydrodynamic model

### Explanations of model

For the one dimensional hydrodynamic models, the forces on netting are calculated based on individual twines. The twines are taken as a cylindrical elements. In practice, the force is usually decomposed into two components: normal drag force ![formula](https://render.githubusercontent.com/render/math?math=\vec{F_n}) and tangential drag force ![formula](https://render.githubusercontent.com/render/math?math=\vec{F_t}) (Cheng *et al*., 2020)

![Fig.5](./figures/Fig.5.png)


![formula](https://render.githubusercontent.com/render/math?math=\quad\Large\vec{F_t}=0.5C_{t}\rho_{w}d_{w}L\left|\vec{U}-\vec{v}\right|^{2}\vec{i_t})


![formula](https://render.githubusercontent.com/render/math?math=\quad\Large\vec{F_n}=0.5C_{n}\rho_{w}d_{w}L\left|\vec{U}-\vec{v}\right|^{2}\vec{i_n})



* ![formula](https://render.githubusercontent.com/render/math?math=\rho_{w}) is the fluid density.

* ![formula](https://render.githubusercontent.com/render/math?math=d_w) is diameter of the cylindrical element.

* ![formula](https://render.githubusercontent.com/render/math?math=L) is the length of the cylindrical element.

* ![formula](https://render.githubusercontent.com/render/math?math=\vec{U}) is the undisturbed incoming flow velocity in the upstream of the net panel.

* ![formula](https://render.githubusercontent.com/render/math?math=\vec{v}) is the velocity of the structure.

* The unit vectors ![formula](https://render.githubusercontent.com/render/math?math=\vec{i_n}) $ and ![formula](https://render.githubusercontent.com/render/math?math=\vec{i_t}) $ which are used to indicate the directions of forces are defined by the following equations (![formula](https://render.githubusercontent.com/render/math?math=\vec{e_i}) is the unit vector of the cylindrical element):


![formula](https://render.githubusercontent.com/render/math?math=\Large\vec{i_t}=\frac{(\vec{U}-\vec{v})\vec{e_i}}{\left|\vec{U}-\vec{v}\right|}\vec{e_i})

![formula](https://render.githubusercontent.com/render/math?math=\Large\vec{i}_{n}=\frac{(\vec{U}-\vec{v})}{\left|\vec{U}-\vec{v}\right|}-\vec{i_t})


* ![formula](https://render.githubusercontent.com/render/math?math=C_n) and ![formula](https://render.githubusercontent.com/render/math?math=C_t) are the normal drag and tangential force coefficients in the one dimension hydrodynamic models, respectively. These force coefficients are usually obtained from experiments that approximate the ideal conditions of a cylinder in an infinite flow field.

### How to implement in code

Example code:
``` python
import hydroModules as hdm

rope = hdm.one_dimensional.line("M4",
                                meshinfo['Lines_rope'],
                                0,
                                caseinfo['Frame']['rope_sec_dia'],
                                caseinfo['Frame']['rope_sec_dia'])
```

> **model_index**: [string] [-]. To indicate the model function, *e.g.*, 'M1', 'M2', 'M3'.
>
> **hydro_element**: [list] [-]. A python list to indicate how the lines are connected, *e.g.*., [[0,1],[1,2],[2,1]]
>
> **solidity**: [float] [-]. The solidity of netting.
>
> **dw0**: [float] [m]. The diameter of the physical net twines. It is used for the hydrodynamic coefficients.
> 
> **dwh**: [float] [m]. The hydrodynamic diameter of the numerical net twines. It is used for the force calculation (reference area).

