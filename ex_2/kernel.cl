__kernel void SAXPY(float a, __global float *X, __global float *Y, ulong size) {
  size_t index = get_global_id(0);
  if (index >= size) return;
  Y[index] += a * X[index];
}
