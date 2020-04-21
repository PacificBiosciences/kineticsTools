from multiprocessing.process import Process, current_process
from kineticsTools.sharedArray import SharedArray

# Test script for making sure that shared memory backed numpy arrays work properly
# FIXME -- migrate to test dir if possible


class Sub(Process):

    def __init__(self, sa):
        Process.__init__(self)
        self.sa = sa
        self.arr = sa.getNumpyWrapper()

    def run(self):
        import time

        print("self.arr[10] = %f, Process = %s" %
              (self.arr[10], current_process()))

        print(self.arr.shape)

        n = self.arr.shape[0] - 1

        print("self.arr[%d] = %f, Process = %s" %
              (n, self.arr[n], current_process()))
        time.sleep(10)


class Test:

    def __init__(self):
        self.sa = SharedArray(dtype="f", shape=50000000)

        self.arr = self.sa.getNumpyWrapper()
        self.arr[:] = 99.0
        self.arr[10] = 100.0

    def start(self):
        proc = Sub(self.sa)
        proc.start()
        proc.join()


if __name__ == "__main__":
    tester = Test()
    tester.start()
