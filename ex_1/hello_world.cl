__kernel void hello_world() {
	int index = get_global_id(0);
	printf("Hello World! My threadId is %d\n", index);
}
