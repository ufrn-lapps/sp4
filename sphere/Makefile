main: sphere reference

sphere: propagation-sphere.c propagation-sphere.h source.c
	gcc source.c propagation-sphere.c -o propagation-sphere -lm
reference: propagation-reference.c propagation-reference.h source.c
	gcc source.c propagation-reference.c -o propagation-reference -lm
clean:
	rm propagation-sphere propagation-reference