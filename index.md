# ScrewShaft as an Elastic body & Coulomb Friction

The theme, which will be completed in about six months, combines NSK's knowledge with threaded shaft elasticity to solve the static force balance, including friction (Coulomb friction only).

## Containts

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=4 orderedList=false} -->


## 1. Discretized representation of beam elements
:

### 1.1. Representation using quaternions
:

#### 1.1.1. Why quaternions?
:

#### 1.1.2. Previous research
:

### 1.2. How to get coordinates in inertial coordinate system
:

### 1.3. Testing the created model
:

#### 1.3.1 Comparison with Eulerian coordinates
:

#### 1.3.2. Representation of a beam in arbitrary coordinates
:

## 2. Discretized representation of a Screw Shaft
:

### 2.1. Not only the nodes, but also the profiles along the way
:

#### 2.1.1. Importance of continuity in dynamic analysis
:

#### 2.1.2. No variable is treated as a constant in the spiral algorithm
:

### 2.2. Mathematical representation of a spiral around a stretched and bent beam
:

### 2.3. Creating an inverse function algorithm
:

### 2.4. Testing the created model
:

#### 2.4.1. Checking the accuracy of the inverse function
:

#### 2.4.2. Is the NSK model correct in the first place? (Confirmation of the Depth of contacting)
:

#### 2.4.3. Is the new model correct? (Confirmation of the Depth of contacting)
:

## 3. The force acting between the ball and the screw shaft including Coulomb friction (in the case of only one ball)
:

### 3.1. What kind of balance is established in the steady state of Coulomb friction?
:

### 3.2. Formulation of nonlinear simultaneous equations
:

## 4. Balance of forces when the screw shaft is treated as an elastic body
:

### 4.1. Modeling of variables in a program
:

#### 4.1.1. Creating a class diagram
:

#### 4.1.2. Creating a test case
:

### 4.2. No Coulomb friction (conventional model)
:

#### 4.2.1. Formulation of nonlinear simultaneous equations
:

#### 4.2.2. Comparison with Bo-Lin model (if possible)
:

### 4.3. With Coulomb friction
:

#### 4.3.1. Formulation of nonlinear simultaneous equations
:




