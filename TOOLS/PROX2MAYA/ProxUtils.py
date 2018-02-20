from math import atan2, asin, pi, sqrt


def quaternion_to_matrix(qs, qx, qy, qz):

    #Make sure we have a unit quaternion
    qnorm = sqrt( qs*qs + qx*qx + qy*qy + qz*qz)

    qs = qs/ qnorm
    qx = qx/ qnorm
    qy = qy/ qnorm
    qz = qz/ qnorm

    R = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    R[0][0] = 1.0 - 2.0 * (qy*qy + qz*qz)
    R[1][1] = 1.0 - 2.0 * (qx*qx + qz*qz)
    R[2][2] = 1.0 - 2.0 * (qy*qy + qx*qx)
    R[1][0] = 2.0 * (qx*qy + qs*qz)
    R[0][1] = 2.0 * (qx*qy - qs*qz)
    R[2][0] = 2.0 * (-qs*qy + qx*qz)
    R[0][2] = 2.0 * (qs*qy + qx*qz)
    R[2][1] = 2.0 * (qz*qy + qs*qx)
    R[1][2] = 2.0 * (qz*qy - qs*qx)
    return R


def euler_angles_from_matrix(R):
    r00 = R[0][0]
    r01 = R[0][1]
    r02 = R[0][2]
    r10 = R[1][0]
    r20 = R[2][0]
    r21 = R[2][1]
    r22 = R[2][2]

    if r20 >= 1.0:
        rz = atan2(-r01, -r02)
        ry = - pi / 2.0
        rx = 0.0
    elif r20 <= -1.0:
        rz = atan2(-r01, r02)
        ry = pi / 2.0
        rx = 0.0
    else:
        rz = atan2(r10, r00)
        ry = asin(-r20)
        rx = atan2(r21, r22)

    # Convert from radians to degrees
    rx = 180.0*rx / pi
    ry = 180.0*ry / pi
    rz = 180.0*rz / pi

    return rx, ry, rz


def make_euler_angles(qs, qx, qy, qz):
    R = quaternion_to_matrix(qs, qx, qy, qz)
    return euler_angles_from_matrix(R)

