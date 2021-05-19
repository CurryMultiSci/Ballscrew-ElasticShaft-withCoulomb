## 1. Discretize of beam elements
Improving Timoshenko's Model Represented by Euler Angles.


### 1.1. Representation using quaternions
Since NSK uses quaternions to handle posture, we will model the quaternion representation in this program as well. Although there is an idea to convert Euler angles and quaternions for each calculation, we will use the direct quaternion method this time because it seems to be computationally expensive and there are already previous studies.

#### 1.1.1. Why quaternions?
Light computational cost. Gimbal lock does not occur.

|                 | Parameters | Memory | Handling | Notes                                                                   |
|-----------------|------------|--------|----------|-------------------------------------------------------------------------|
| Quaternion      | 4          | 0.5    | 1        | Has a good balance of abilities                                         |
| Eular           | 3          | 1      | 0        | The gradient per unit attitude change differs depending on the attitude |
| Rotation Matrix | 9          | 0      | 0.5      |                                                                         |


```python
import numpy as np
# from numpy import quaternion

```


```python


```

#### 1.1.2. Previous research

Two types of modeling of beams by quaternions have been presented by the same laboratory.

|Kinematically Exact Beam Finite Elements Based on Quaternion Algebra (2014)|A consistent strain-based beam element with quaternion representation of rotations (2020)|
|---|---|
|![1‗2014](1_2014.png)|![1_2020](1_2020.png)|
|Eva Zupan, Miran Saje, **Dejan Zupan**|Damjan Loli, **Dejan Zupan**, Miha Brojan|
|University of Ljubljana|University of Ljubljana|
|Simple|Complex|
|displacement-based|strain-based|
|12 DOFs per element|26 DOFs per element|
|**to small shear deformations**|to larger shear deformations|

In order to deal with the deformation of the ballscrew this time, we do not need to consider "larger shear deformations", so I will deal with the **2014 model**.



```python

```
