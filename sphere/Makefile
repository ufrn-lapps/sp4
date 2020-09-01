main: sphere reference sphere-3d reference-3d

# 2D Version
sphere: 2d/propagation-sphere.c 2d/propagation-sphere.h source.c
	gcc source.c 2d/propagation-sphere.c -o 2d/propagation-sphere -lm
reference: 2d/propagation-reference.c 2d/propagation-reference.h source.c
	gcc source.c 2d/propagation-reference.c -o 2d/propagation-reference -lm

# 3D Version
reference-3d: 3d/propagation-reference-3d.c 3d/propagation-reference-3d.h source.c
	gcc -Wall source.c 3d/propagation-reference-3d.c -o 3d/propagation-reference-3d -lm

sphere-3d: 3d/propagation-sphere-3d.c 3d/propagation-sphere-3d.h source.c
	gcc source.c 3d/propagation-sphere-3d.c -o 3d/propagation-sphere-3d -lm

clean:
	rm 2d/propagation-sphere 2d/propagation-reference 3d/propagation-reference-3d 3d/propagation-sphere-3d