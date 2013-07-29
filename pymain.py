from ctypes import CDLL, c_double

import random
import math
import threading


def main():
    dll = CDLL('smallpt.so')
    dll.addSphere.argtypes = [c_double, c_double, c_double, c_double]
    dll.addLight.argtypes = [c_double, c_double, c_double, c_double]

    for i in range(100):
        cx = cy = cz = 50
        r = random.random() * 50
        theta = random.random() * 2 * math.pi
        phi = random.random() * math.pi

        params = (
            random.random() * 5,
            cx + r * math.sin(phi) * math.cos(theta),
            cy + r * math.sin(phi) * math.sin(theta),
            cz + r * math.cos(phi))
        if random.random() < 0.1:
            dll.addLight(*params)
        else:
            dll.addSphere(*params)
    t = threading.Thread(target=lambda: dll.render("thing.ppm", 1024, 768, 32))
    t.daemon = True
    t.start()
    while t.is_alive():
        # running in a separate daemon thread with a timeout on the join
        # allows Ctrl-C interruption
        t.join(timeout=10)


if __name__ == '__main__':
    main()
