# TEST SCRIPT FOR DNA_SEQUENCER

# Difference()
echo "d\nfile\ntest\nd\nq\n" >  input.txt
make irun

echo "d\nfile\ntest\nm\n100\n50\n3\n2\n10\nq\n" >  input.txt
make irun

echo "d\nstdout\nd\nq\n" >  input.txt
make irun

# Test()
echo "t\n10\n\nd\nq\n" >  input.txt
make irun
