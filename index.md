# ScrewShaft as an Elastic body & Coulomb Friction

The theme, which will be completed in about six months, combines NSK's knowledge with threaded shaft elasticity to solve the static force balance, including friction (Coulomb friction only).

## Containts

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=4 orderedList=false} -->


## Discretized representation of beam elements
:

### Representation using quaternions
:

#### Why quaternions?
:

#### Previous research
:

### How to get coordinates in inertial coordinate system
:

### Testing the created model
:

#### Comparison with Eulerian coordinates
:

#### Representation of a beam in arbitrary coordinates
:

## Discretized representation of a Screw Shaft
:

### Not only the nodes, but also the profiles along the way
:

#### Importance of continuity in dynamic analysis
:

#### No variable is treated as a constant in the spiral algorithm
:

### Mathematical representation of a spiral around a stretched and bent beam
:

### Creating an inverse function algorithm
:

### Testing the created model
:

#### Checking the accuracy of the inverse function
:

#### Is the NSK model correct in the first place? (Confirmation of the Depth of contacting)
:

#### Is the new model correct? (Confirmation of the Depth of contacting)

:

## The force acting between the ball and the screw shaft including Coulomb friction (in the case of only one ball)
:

### What kind of balance is established in the steady state of Coulomb friction?
:

### Formulation of nonlinear simultaneous equations
:

## Balance of forces when the screw shaft is treated as an elastic body
:

### Modeling of variables in a program
:

#### Creating a class diagram
:

#### Creating a test case
:

### No Coulomb friction (conventional model)
:

#### Formulation of nonlinear simultaneous equations
:

#### Comparison with Bo-Lin model (if possible)
:

### With Coulomb friction
:

#### Formulation of nonlinear simultaneous equations
:




