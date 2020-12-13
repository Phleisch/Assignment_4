__kernel void hello_world() {
	int index = get_global_id(0) + 16*get_global_id(1) + 64*get_global_id(2);
	printf("Hello World! My threadId is %d\n", index);
}
