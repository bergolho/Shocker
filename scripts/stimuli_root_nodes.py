import sys
import numpy as np

def main():
	
	if len(sys.argv) != 2:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file> = Input file with the root node locations")
		print("-------------------------------------------------------------------------")
		return 1

	filename = sys.argv[1]
	data = np.genfromtxt(filename)
	for i in range(len(data)):
		x, y, z, lat = data[i][0], data[i][1], data[i][2], data[i][3]
		print("[stim_root_node_%d]" % (i))
		print("start = %.1lf" % (lat))
		print("duration = 4.0")
		print("current = -53.0")
		print("center_x = %g" % (x))
		print("center_y = %g" % (y))
		print("center_z = %g" % (z))
		print("radius = 2000.0")
		print("main_function = stim_sphere")
		print()

if __name__ == "__main__":
	main()
