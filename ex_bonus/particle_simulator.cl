typedef struct {
  float3 position;
  float3 velocity;
} Particle;

__kernel void particle_simulator(__global Particle *particles, float3 field, int num_particles, int num_iterations) {
	int threadId = 0;
	if (threadId > num_particles) return;
	Particle particle = particles[threadId];
	for (int i = 0; i < num_iterations; i++) {
		particle.velocity.x = particle.velocity.x + field.x;
		particle.velocity.y = particle.velocity.y + field.y;
		particle.velocity.z = particle.velocity.z + field.z,
		
		particle.position.x = particle.position.x + particle.velocity.x;
		particle.position.y = particle.position.y + particle.velocity.y;
		particle.position.z = particle.position.z + particle.velocity.z;
	}
	particles[threadId] = particle;
}

