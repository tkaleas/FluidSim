Smoke Simulation

I basically completed the assignment as directed.

In order to add solid cells in the program I added a new GridData mS to store whether or not there is a solid in each cell.
I also colorized the cells in the renderer so you can visualize the solid cells, which solid grey cubes.

However I added a preconditioner and a modified cg_solve to MACGrid.h/.cpp so that I could apply it and get a speedup.

The preconditioner I initially used was the Incompletely Conjugate Gradient, or IC(0), which got me a speedup of 1/2 of the iterations needed.
However I modified it again and added the MIC(0), which then got me a speed up which was a little bit less than half of the iterations needed (however in this
case it scales up so that I get a better speedup for systems that are very large and have a lot of grid cells and have a lot going on in them (lots of vortices and a lot of movement inbetween cells).



Answers to Theoretical Questions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
Soft Constraints
----------------

1) The first approach, which we used in our project, was a soft constraint,
in switch when the particle went "off" track or tried to penetrate an object, we 
applied a force on it to try and place it back on track from where it was
leaving the constrained path. We basically used collision detection to determine
whether or not to apply a constraining force on the particle.

This approach is simple and easy to implement, and has a low computational cost, but the particle will most likely
not keep directly at the constraint, because the forces are acting on changes in the particle's
position rather than the positions themselves. Thus it is impossible to achieve accuracy without applying
and maintaining stiff spring forces on the particles.

2) 
A similar approach using soft constraints is to create "behavior funcions" whose equations determine where a particle should 
or should not be. We define a scalar potential energy function, and then determine the forces on the particle so that the behavior function
meets some constraints (usually that it is equal to 0), and update every time step. This is a little bit better than the method we used in our
strict penalty method, and in a sense we used a type of it, however as a soft constraint the method still suffers from the fact that the
constraint will never be full and will still be innacurate.

Hard Constraints
-----------------
3) Another method is representing the constraint using an implicit equation for the constraint,
which we can then use to calculate the constraining force we should apply on the particle
to make sure the constraint is met. 

This method has the advantage of making sure that our constraints are met in the system with
none of the drift and innaccuracy of the soft constraints,
and through some linear algebra we can solve a large system of constraints fairly efficiently.

4) Another similar hard constraint method would be to use parametric equations to reduce the "degrees of freedom" 
and as well as the things we have to solve for in order to calculate the constraining forces. This is helpful because it means
that we have fewer things so solve for, as well as the fact the constraints are always met. However it is difficult to combine
multiple constraints as well as formulate the constraints using parametric constraints.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
2) One of the big problems that we might have run into with working with a hard constraint system is accounting for a dynamic 
and changing constraint system as well as a changing particle system. It is extremely difficult to try and hand model a large set of constraints
for instance in a system like ours. What we would need is some way to create modular constraints that we can combine together and can
model "on the fly" instead of having only 1. A solution to this would be to create a fragmented system, with indiviual contraints similar to the way we create individual forces,
 where each constraint evaluates only itself and its derivatives, and then all of the constraints are combined later on.
 
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 3) Repulsion forces are more efficient than constraining dynamics because within the cloth, we will have many many penetrations that are constantly occuring between the cloth and itself in our case,
 so it is a lot easier computationally if we just apply soft repulsion forces on the particles rather than solve for constraining forces, which may bring other problems. For instance, if we get stuck in a 
time step which somehow a lot of interpenetration can occur, it might mean that constraining forces may keep the particle embedded in another part of the cloth. However with penalty forces for cloth/cloth penetration
the spring forces that we will apply will change the courses of both cloth particles to ensure that they do not interpenetrate.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
4)  As far as explicit integration goes, there are big problems with the stability of the simulation. Meaning that if your step size is too big, the simulation could blow up as it runs wild and misses some asymptoptic behavior (perhaps)
 Reducing the step size could remedy the situation however when you add things like constraining dynamics the system could become unstable regardless of how small the step size is.
 
 5) The midpoint method is within O(h^3)  error where h is the time step.
 The derivation of this is from the Taylor Series:
 
 x(t0+h)=x(t0) + hx(t0) + h^2x(t0)/2! + h^3x(t0)/3! ...

so if we drop this last noted term we will end up with an error that is equal to some h^3 * x(t0)/3! error away from the actual value, or O(h^3). 
 
 6) Explicit integration does not work well for cloth systems, because in cloth we use a lot of stiff springs to simulate the cloth, and explicit integration doesn't work well for stiff systems.
 A stiff system is one which the function changes a great deal on a short interval in some cases, or in some cases it changes extremely little. In our case for cloth, it means that we have some springs which do not like to change,
 and other springs which apply extremely strong forces in some directions in order to keep 2 particles with springs between them at a certain distance, So we can have large changes in forces over a very short timestep 
 Explicit integration systems cannot handle stiff systems because given a certain time step, they may skip some important part of a function and make an extreme large forces instead of asympoptically arriving at a smaller one, which
 leads to the system becoming unstable and blowing up.
 
 7) In the method Baraff describes in his notes, he already needs to solve a linear system at each step. If we moved this to a higher order integration, we may need to solve the system for a magnitude higher so it may dramatically increase
 the computational cost of integration. However, this system, as he describes, will most likely be sparse as well as the fact that we gain better accuracy by using a higher order term.
 
 8) The main contribution of "Large Steps in Cloth Simulation" is using a variety of techniques to make it so that you can use large time steps in order to model cloth. The main part of this way to use larger timesteps to model cloth was
 using IMPLICIT INTEGRATION.
 
 9) The main contribution of "Robust Treatment of Collisions..." was that it described methods to handle all kinds of collision detection and response in cloth simulation, and it was the first paper to resolve all of the collisions in cloth.
 (Plus the videos that came out of it looked really really good).