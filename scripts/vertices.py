import sys

def main():
	
	if len(sys.argv) != 2:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <N>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<N> = Number of vertices to plot")
		print("-------------------------------------------------------------------------")
		return 1

	N = int(sys.argv[1])
	print("VERTICES %u %u" % (N,N*2))
	for i in range(N):
		print("1 %u" % (i))


if __name__ == "__main__":
	main()
