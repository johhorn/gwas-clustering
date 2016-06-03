from subprocess import Popen, PIPE
import sys, getopt, os

def main(argv):
   	folder = None
   	perlfile = None
   	rsquare = '0'
   	error = sys.argv[0] + ' -p <perlfile> -f <folder> -r <rsquare>'

   	try:
		opts, args = getopt.getopt(argv,"hp:f:r:",["pfile=", "folder=", "rsquare="])
	except getopt.GetoptError:
		print error
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print error
			sys.exit()		
		elif opt in ("-p", "--pfile"):
			perlfile = arg
		elif opt in ("-f", "--folder"):
			folder = arg
		elif opt in ("-r", "--rsquare"):
			rsquare = arg

	if not folder or not perlfile:
		print error
		sys.exit()

	path = folder
	inputfiles = []
	outputfiles = []

	for root, dirs, files in os.walk(path, topdown=False):
		for name in files:

			if name.startswith('out') or name.startswith('long'):
				continue

			inputfiles.append(os.path.join(root, name))
			out = "out_" + rsquare + '_' + name
			outputfiles.append(os.path.join(root, out))

	for i in range(0, len(inputfiles)):
		inputfile = inputfiles[i]
		outputfile = outputfiles[i]
		
		print 'inputfile', inputfile
		print 'outputfile', outputfile

		pipe = Popen(["perl", perlfile, inputfile, outputfile, rsquare], stdin=PIPE, stdout=PIPE)
		result = pipe.stdout.read()
		print result

if __name__ == "__main__":
	main(sys.argv[1:])