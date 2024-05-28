#include "acclaim/motion.h"
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include "simulation/kinematics.h"
#include "util/helper.h"

namespace acclaim {
Motion::Motion() {}
Motion::Motion(const util::fs::path &amc_file, std::unique_ptr<Skeleton> &&_skeleton) noexcept
    : skeleton(std::move(_skeleton)) {
    postures.reserve(1024);
    if (!this->readAMCFile(amc_file)) {
        std::cerr << "Error in reading AMC file, this object is not initialized!" << std::endl;
        std::cerr << "You can call readAMCFile() to initialize again" << std::endl;
        postures.resize(0);
    }
}

const std::unique_ptr<Skeleton> &Motion::getSkeleton() const { return skeleton; }

Motion::Motion(const Motion &other) noexcept
    : skeleton(std::make_unique<Skeleton>(*other.skeleton)), postures(other.postures) {}

Motion::Motion(Motion &&other) noexcept : skeleton(std::move(other.skeleton)), postures(std::move(other.postures)) {}

Motion::Motion(const Motion &other, int startIdx, int endIdx) : skeleton(std::make_unique<Skeleton>(*other.skeleton)) {
    if (startIdx < 0 || startIdx >= other.postures.size() || endIdx <= startIdx || endIdx > other.postures.size()) {
        throw std::out_of_range("Invalid start or end index");
    }

    postures.assign(other.postures.begin() + startIdx, other.postures.begin() + endIdx);
}

Motion &Motion::operator=(const Motion &other) noexcept {
    if (this != &other) {
        skeleton.reset();
        skeleton = std::make_unique<Skeleton>(*other.skeleton);
        postures = other.postures;
    }
    return *this;
}

Motion &Motion::operator=(Motion &&other) noexcept {
    if (this != &other) {
        skeleton = std::move(other.skeleton);
        postures = std::move(other.postures);
    }
    return *this;
}

void Motion::remove(int begin, int end) {
    assert(end >= begin);
    postures.erase(postures.begin() + begin, postures.begin() + end);
}

void Motion::concatenate(Motion &m2) {
    for (int i = 0; i < m2.getFrameNum(); i++)
        postures.insert(postures.end(), m2.postures[i]);
}

Posture& Motion::getPosture(int FrameNum) { 
    return postures[FrameNum]; 
}

std::vector<Posture> Motion::getPostures() { 
    return postures; 
}

void Motion::setPosture(int FrameNum, const Posture &InPosture) {
    postures[FrameNum] = InPosture; 
}

int Motion::getFrameNum() const { return static_cast<int>(postures.size()); }

void Motion::forwardkinematics(int frame_idx) {
    kinematics::forwardSolver(postures[frame_idx], skeleton->getBonePointer(0));
    skeleton->setModelMatrices();
}

void Motion::transform(Eigen::Vector4d& newFacing, const Eigen::Vector4d &newPosition) {
    // *TODO*
    // Task: Transform the whole motion segment so that the root bone of the first posture(first frame)
    //       of the motion is located at newPosition, and its facing be newFacing.
    //       The whole motion segment must remain continuous.
    
    if (postures.empty()) return;

    Eigen::Vector4d initPos = postures[0].bone_translations[0];
    Eigen::Vector4d delta = newPosition - initPos;
    
    Eigen::Matrix3d newRotMat = util::rotateDegreeZYX(newFacing).normalized().toRotationMatrix();
    Eigen::Matrix3d initRotMat = util::rotateDegreeZYX(postures[0].bone_rotations[0]).normalized().toRotationMatrix();

    double thetaY = atan2(newRotMat(2, 2), newRotMat(0, 2)) - atan2(initRotMat(2, 2), initRotMat(0, 2));

    Eigen::Matrix3d rotMatR = Eigen::AngleAxisd(thetaY, Eigen::Vector3d::UnitY()).toRotationMatrix();

    for (auto &posture : postures) {
        Eigen::Vector3d translate = rotMatR *(posture.bone_translations[0].head<3>() - initPos.head<3>()) + initPos.head<3>() +
            delta.head<3>();
        posture.bone_translations[0] = Eigen::Vector4d(translate[0], translate[1], translate[2], 0.0);
            
        
        Eigen::Matrix3d currentRotationMatrix = util::rotateDegreeZYX(posture.bone_rotations[0]).normalized().toRotationMatrix();
        Eigen::Matrix3d newRotationMatrix = rotMatR * currentRotationMatrix;
        //Eigen::Vector3d newRotationEuler = util::Quater2EulerAngle(Eigen::Quaterniond(newRotationMatrix));
        Eigen::Vector3d newRotationEuler = newRotationMatrix.eulerAngles(0, 1, 2);
        posture.bone_rotations[0] =
            util::toDegree(Eigen::Vector4d(newRotationEuler[0], newRotationEuler[1], newRotationEuler[2], 0.0));
    }
}
    
Motion blend(Motion bm1, Motion bm2, const std::vector<double> &weight) {
    // *TODO*
    // Task: Return a motion segment that blends bm1 and bm2.
    //       bm1: tail of m1, bm2: head of m2
    //       You can assume that m2's root position and orientation is aleady aligned with m1 before blending.
    //       In other words, m2.transform(...) will be called before m1.blending(m2, blendWeight, blendWindowSize) is
    //       called

    int blendWindowSize = bm1.getFrameNum();
    int boneNum = bm1.getPosture(0).bone_translations.size();
    Motion blendedMotion = bm2;

    for (int i = 0; i < blendWindowSize; i++) {
        Posture resultPosture = blendedMotion.getPosture(i);

        std::cout << weight[i] << '\n';
        resultPosture.bone_translations[0] = bm1.getPosture(i).bone_translations[0] * weight[i] +
                                             bm2.getPosture(i).bone_translations[0] * (1.0 - weight[i]);
        // resultPosture.bone_rotations[0]    = bm1.getPosture(i).bone_rotations[0] * weight[i] + 
        //                                      bm2.getPosture(i).bone_rotations[0] * (1.0 - weight[i]);

        // Interpolate rotations (using quaternion slerp)
        Eigen::Quaterniond q1(bm1.getPosture(i).bone_rotations[0]);
        Eigen::Quaterniond q2(bm2.getPosture(i).bone_rotations[0]);
        Eigen::Quaterniond qInterpolated = q1.slerp(weight[i], q2);

        // Convert interpolated quaternion back to euler angles
        Eigen::Vector3d newRotationEuler = util::Quater2EulerAngle(qInterpolated);
        resultPosture.bone_rotations[0].head<3>() = newRotationEuler;

        blendedMotion.setPosture(i, resultPosture);
    }

    return blendedMotion;
}

Motion Motion::blending(Motion &m2, const std::vector<double> &blendWeight, int blendWindowSize) {
    // do blending
    Motion bm1(*this, this->getFrameNum() - blendWindowSize, this->getFrameNum());
    Motion bm2(m2, 0, blendWindowSize);
    Motion bm = blend(bm1, bm2, blendWeight);
    return bm;
}


bool Motion::readAMCFile(const util::fs::path &file_name) {
    // Open AMC file
    std::ifstream input_stream(file_name);
    // Check if file successfully opened
    if (!input_stream) {
        std::cerr << "Failed to open " << file_name << std::endl;
        return false;
    }
    // There are (NUM_BONES_IN_ASF_FILE - 2) moving bones and 2 dummy bones (lhipjoint and rhipjoint)
    int movable_bones = skeleton->getMovableBoneNum();
    // Ignore header
    input_stream.ignore(1024, '\n');
    input_stream.ignore(1024, '\n');
    input_stream.ignore(1024, '\n');
    int frame_num;
    std::string bone_name;
    while (input_stream >> frame_num) {
        auto &&current_posture = postures.emplace_back(skeleton->getBoneNum());
        for (int i = 0; i < movable_bones; ++i) {
            input_stream >> bone_name;
            const Bone &bone = *skeleton->getBonePointer(bone_name);
            int bone_idx = bone.idx;
            Eigen::Vector4d bone_rotation = Eigen::Vector4d::Zero();
            Eigen::Vector4d bone_translation = Eigen::Vector4d::Zero();
            if (bone.doftx) {
                input_stream >> bone_translation[0];
            }
            if (bone.dofty) {
                input_stream >> bone_translation[1];
            }
            if (bone.doftz) {
                input_stream >> bone_translation[2];
            }
            if (bone.dofrx) {
                input_stream >> bone_rotation[0];
            }
            if (bone.dofry) {
                input_stream >> bone_rotation[1];
            }
            if (bone.dofrz) {
                input_stream >> bone_rotation[2];
            }
            current_posture.bone_rotations[bone_idx] = std::move(bone_rotation);
            current_posture.bone_translations[bone_idx] = std::move(bone_translation);
            if (bone_idx == 0) {
                current_posture.bone_translations[bone_idx] *= skeleton->getScale();
            }
        }
    }
    input_stream.close();
    std::cout << frame_num << " samples in " << file_name.string() << " are read" << std::endl;
    return true;
}
void Motion::render(graphics::Program *program) const { skeleton->render(program); }

}  // namespace acclaim
