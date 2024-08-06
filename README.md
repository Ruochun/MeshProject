 To use it:
 g++ MeshProject.cpp -o test


 Discussion:
 1. The geometry hierarchy usually requires point, edge and faces to be stored. I have them stored in vectors for fast access. 
 2. First we'll need to check the Jacobian of mesh/aspect ratio of the elements after the operation, and this is concerning if such operations happened near the shape edges. Then we should ensure that the topology of the mesh is still valid.