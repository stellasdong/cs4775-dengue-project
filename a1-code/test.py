score_mtx = {"A": {"A": 91, "C": -114, "T": -123, "G": -31}, 
	"C": {"A": -114, "C": 100, "T": -31, "G": -125}, 
	"T": {"A": -123, "C": -31, "T": 91, "G": -114}, 
	"G": {"A": -31, "C": -125, "T": -114, "G": 100}}


# Actual result (Please substitute these to your own alignment)
a_x = "GGGTGGGAAAATAGACCAATAGGCAGAGAGAGTCAGTGCCTATCAGAAACCCAAGAGTCTTCTCTGTCTCCACATGCCCAGTTTCTATTGGTCTCCTTAAACCTGTCTTGTAACCTTGATA"
a_y = "AAAGGGAAACATAGA-CAGGGGACACTCAAAGTTAGTGCCTGCTGGAAAGCAGA------CCTCTGTCTCCAAGCACCCAACTTCTACTTGT-------GAGCTGCCTTGTAACCTGGATA"


d = 430
e = 30

current_status = "match_mismatch"
total_score = 0
for i in range(len(a_x)):
	if current_status == "match_mismatch":
		if a_x[i]=="-":
			total_score = total_score - d
			current_status = "gap_x"
		elif a_y[i]=="-":
			total_score = total_score - d
			current_status = "gap_y"
		else:
			total_score = total_score + score_mtx[a_x[i]][a_y[i]]
	elif current_status == "gap_x":
		if a_x[i]=="-":
			total_score = total_score - e
		else:
			total_score = total_score + score_mtx[a_x[i]][a_y[i]]
			current_status = "match_mismatch"
	elif current_status == "gap_y":
		if a_y[i]=="-":
			total_score = total_score - e
		else:
			total_score = total_score + score_mtx[a_x[i]][a_y[i]]
			current_status = "match_mismatch"
	else:
		print("?")



print(total_score)