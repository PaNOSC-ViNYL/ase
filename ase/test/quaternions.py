import numpy as np
from ase.quaternions import Quaternion


def axang_rotm(u, theta):

    u = np.array(u).astype(np.float)
    u /= np.linalg.norm(u)

    # Cross product matrix for u
    ucpm = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])

    # Rotation matrix
    rotm = (np.cos(theta) * np.identity(3) + np.sin(theta) * ucpm +
            (1 - np.cos(theta)) * np.kron(u[:, None], u[None, :]))

    return rotm


def rand_rotm():
    """Axis & angle rotations."""
    u = np.random.random(3)
    theta = np.random.random() * np.pi * 2

    return axang_rotm(u, theta)


def eulang_rotm(a, b, c, mode='zyz'):

    rota = axang_rotm([0, 0, 1], a)
    rotc = axang_rotm([0, 0, 1], c)

    if mode == 'zyz':
        rotb = axang_rotm([0, 1, 0], b)
    elif mode == 'zxz':
        rotb = axang_rotm([1, 0, 0], b)

    return np.dot(rotc, np.dot(rotb, rota))

# First: test that rotations DO work
for i in range(10):
    # 10 random tests

    rotm = rand_rotm()

    q = Quaternion.from_matrix(rotm)

    # Now test this with a vector
    v = np.random.random(3)

    vrotM = np.dot(rotm, v)
    vrotQ = q.rotate(v)

    assert np.allclose(vrotM, vrotQ)

# Second: test the special case of a PI rotation

rotm = np.identity(3)
rotm[:2, :2] *= -1               # Rotate PI around z axis

q = Quaternion.from_matrix(rotm)

assert not np.isnan(q.q).any()

# Third: test compound rotations and operator overload
for i in range(10):

    rotm1 = rand_rotm()
    rotm2 = rand_rotm()

    q1 = Quaternion.from_matrix(rotm1)
    q2 = Quaternion.from_matrix(rotm2)

    # Now test this with a vector
    v = np.random.random(3)

    vrotM = np.dot(rotm2, np.dot(rotm1, v))
    vrotQ = (q2 * q1).rotate(v)

    assert np.allclose(vrotM, vrotQ)

# Fourth: test Euler angles
for mode in ['zyz', 'zxz']:
    for i in range(10):

        abc = np.random.random(3)*2*np.pi
        v2 = np.random.random((2, 3))  # Two random vectors to rotate rigidly

        q_eul = Quaternion.from_euler_angles(*abc, mode=mode)
        rot_eul = eulang_rotm(*abc, mode=mode)

        v2_q = np.array([q_eul.rotate(v) for v in v2])
        v2_m = np.array([np.dot(rot_eul, v) for v in v2])

        assert np.allclose(v2_q, v2_m)
