#include "integrator.h"

#include "configs.h"

void ExplicitEuler::integrate(const std::vector<Particles *> &particles, std::function<void(void)>) const {
  // TODO: Integrate velocity and acceleration
  //   1. Integrate velocity.
  //   2. Integrate acceleration.
  //   3. You should not compute position using acceleration. Since some part only update velocity. (e.g. impulse)
  // Note:
  //   1. You don't need the simulation function in explicit euler.
  //   2. You should do this first because it is very simple. Then you can chech your collision is correct or not.
  //   3. This can be done in 5 lines. (Hint: You can add / multiply all particles at once since it is a large matrix.)

  particles[0]->position() += deltaTime * particles[0]->velocity();
  particles[0]->velocity() += deltaTime * particles[0]->acceleration();
}

void ImplicitEuler::integrate(const std::vector<Particles *> &particles,
                              std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Integrate velocity and acceleration using explicit euler to get Xn+1.
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.

  Eigen::Matrix4Xf xn = particles[0]->position();
  Eigen::Matrix4Xf vn = particles[0]->velocity();

  particles[0]->position() += deltaTime * particles[0]->velocity();
  particles[0]->velocity() += deltaTime * particles[0]->acceleration();

  simulateOneStep();

  particles[0]->position() = xn + deltaTime * particles[0]->velocity();
  particles[0]->velocity() = vn + deltaTime * particles[0]->acceleration();
}

void MidpointEuler::integrate(const std::vector<Particles *> &particles,
                              std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Integrate velocity and acceleration using explicit euler to get Xn+1.
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.

  Eigen::Matrix4Xf xn = particles[0]->position();
  Eigen::Matrix4Xf vn = particles[0]->velocity();

  particles[0]->position() += deltaTime * particles[0]->velocity() / 2;
  particles[0]->velocity() += deltaTime * particles[0]->acceleration() / 2;

  simulateOneStep();

  particles[0]->position() = xn + deltaTime * particles[0]->velocity();
  particles[0]->velocity() = vn + deltaTime * particles[0]->acceleration();
}

void RungeKuttaFourth::integrate(const std::vector<Particles *> &particles,
                                 std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Compute k1, k2, k3, k4
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.

  Eigen::Matrix4Xf xn = particles[0]->position();
  Eigen::Matrix4Xf vn = particles[0]->velocity();

  Eigen::Matrix4Xf k1_pos = deltaTime * particles[0]->velocity();
  Eigen::Matrix4Xf k1_vel = deltaTime * particles[0]->acceleration();
  particles[0]->position() = xn + k1_pos / 2;
  particles[0]->velocity() = vn + k1_vel / 2;
  simulateOneStep();

  Eigen::Matrix4Xf k2_pos = deltaTime * particles[0]->velocity();
  Eigen::Matrix4Xf k2_vel = deltaTime * particles[0]->acceleration();
  particles[0]->position() = xn + k2_pos / 2;
  particles[0]->velocity() = vn + k2_vel / 2;
  simulateOneStep();

  Eigen::Matrix4Xf k3_pos = deltaTime * particles[0]->velocity();
  Eigen::Matrix4Xf k3_vel = deltaTime * particles[0]->acceleration();
  particles[0]->position() = xn + k3_pos;
  particles[0]->velocity() = vn + k3_vel;
  simulateOneStep();

  Eigen::Matrix4Xf k4_pos = deltaTime * particles[0]->velocity();
  Eigen::Matrix4Xf k4_vel = deltaTime * particles[0]->acceleration();

  particles[0]->position() = xn + (k1_pos + 2 * k2_pos + 2 * k3_pos + k4_pos) / 6;
  particles[0]->velocity() = vn + (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel) / 6;
}
