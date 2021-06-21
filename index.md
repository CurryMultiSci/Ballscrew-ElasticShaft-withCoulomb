# Screw Shaft as an Elastic body & Coulomb Friction

The theme, which will be completed in about six months, combines NSK's knowledge with threaded shaft elasticity to solve the static force balance, including friction (Coulomb friction only).

![intro Image 0](/_0_intro/overview.svg)


# Containts

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=4 orderedList=false} -->



## 1. Discretize of beam elements
Improving Timoshenko's Model Represented by Euler Angles.

![1-1](/_1_Discretize_of_beam_elements/1_2014.png)

[Show more](_1_Discretize_of_beam_elements/index.md)

### 1.1. Representation using quaternions
Since NSK uses quaternions to handle posture, we will model the quaternion representation in this program as well. Although there is an idea to convert Euler angles and quaternions for each calculation, we will use the direct quaternion method this time because it seems to be computationally expensive and there are already previous studies.

#### 1.1.1. Why quaternions?
Light computational cost. Gimbal lock does not occur.

#### 1.1.2. Previous research
Lolić, Damjan, Dejan Zupan, and Miha Brojan. "A consistent strain-based beam element with quaternion representation of rotations." Computational Mechanics (2020): 1-16.

### 1.2. How to get coordinates in inertial coordinate system
Since a piece of the screw shaft can move to any position in the inertial coordinate system, its node must correspond to an arbitrary position.

### 1.3. Testing the created model
  

#### 1.3.1 Comparison with Eulerian coordinates
Comparison with previous model.

#### 1.3.2. Representation of a beam in arbitrary coordinates
Place the beam element at an arbitrary position on the inertial coordinate system and check if the calculation does not go wrong.


## 2. Discretize of Screw Shaft
  
![Image2](/_2_Discretize_of_Screw_Shaft/Picture2.svg)


### 2.1. Not only the coordinates of the nodes, but also the profiles along the way
In this case, not only the nodes at both ends but also the center profile is important. There are two reasons for this. One is that the dynamic analysis needs to be continuous, so not only the coordinates of the nodes at both ends, but also the center profile is important. The other reason is that the NSK algorithm does not treat all variables (such as the spiral phase angle) as constants, so it is not possible to locate the nodes in the first place.


### 2.2. Mathematical representation of a spiral around a stretched and bent beam
First, we need to determine the equation of the spiral surrounding the stretched and bent beam, and the Frenet coordinate system.


### 2.3. Creating an inverse function algorithm
Once the Frenet coordinates are determined, the inverse function algorithm is needed. We will have to start from scratch again, but we expect to be able to do so by using quadratic approximation.

### 2.4. Testing the created model
  

#### 2.4.1. Checking the accuracy of the inverse function
Perform a forward transformation → inverse transformation to check if the point returns to its position in the original inertial coordinate system.


#### 2.4.2. Is the NSK model correct in the first place? (Confirmation of the Depth of contacting)
Actually, this part was not fully verified. The position on the right angle section of the groove is calculated correctly, but it is not clear if the contact point position rides on that plane. First of all, let's check this using a simple model.


#### 2.4.3. Is the new model correct? (Confirmation of the Depth of contacting)
Similar to the above method, we will now verify the model with a stretched and bent beam model.


## 3. Friction between ball and spiral
With the current progress of NSK, friction can only be taken into account in dynamic analysis, so static load balancing cannot be obtained and is not suitable to be combined with the current beam model. In addition, users in the company have complained about the calculation speed and convergence judgment. Therefore, the first step is to formulate the balance including the Coulomb friction.

![3-1](/_3_Friction_between_ball_and_spiral/Picture5.svg)

### 3.1. Contact with Spiral Torus

[Show more](_3_Friction_between_ball_and_spiral/_3_1_Contact_with_Spiral_Torus/index.html)


### 3.2. What kind of balance is established in the steady state of Coulomb friction?
When the ball is in a steady state, the force and torque of the ball should be zero (strictly speaking, there are centrifugal force and gyroscopic moment, but their effects are small and should be ignored), but it is very difficult to find the force balance because of the complexity of the spiral shape. By looking at the state of each contact point in the dynamic analysis in detail, we can understand what is happening.

### 3.3. Formulation of nonlinear simultaneous equations
If the mechanism of Coulomb friction force generation is known, it can be formulated, so a simultaneous equation can be formulated and solved by a nonlinear solver. To check if the calculation is correct, the steady state of the dynamic analysis and the position are matched.


## 4. Balance of forces when the screw shaft is treated as an elastic body
  

### 4.1. Modeling of variables in a program
In order to cut the screw shaft as a piece, the contact point position is necessary. However, in this model, the contact point position on the screw side is treated as a variable. Therefore, it is expected to be very difficult to model (so far I have not been able to imagine how to do it). If the balls are completely free, that will not uniquely determine the solution, so I plan to fix one degree of freedom for the relative position as seen from the nut and give two degrees of freedom.


#### 4.1.1. Express the discretization of the screw axis in terms of variables
  

#### 4.1.2. Creating a test case
  

### 4.2. No Coulomb friction (conventional model)
  

#### 4.2.1. Formulation of nonlinear simultaneous equations
  

#### 4.2.2. Comparison with Bo-Lin model (if possible)
  

### 4.3. With Coulomb friction
  

#### 4.3.1. Formulation of nonlinear simultaneous equations
  




