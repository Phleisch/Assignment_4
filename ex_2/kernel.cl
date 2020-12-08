__kernel void SAXPY(float a, __global float *X, __global float *Y) {      
  int index = get_global_id(0);		             
  if (index >= 10000) return;         
  Y[index] += a * X[index];                     
}														                     
