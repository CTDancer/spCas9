import os
import glob
import pandas as pd

number_list = []
score_list = []
with open('/home/ubuntu/scratch/tongchen/spCas9/good_results.txt', 'r') as file:
    lines = file.readlines()

for line in lines:
    number, score = line.split(' | ')
    number_list.append(number)
    score_list.append(score)

print(number_list)
assert len(number_list) == len(score_list)

sequence_list = []
length_list = []
for num in number_list:
    csv_file_path = f'/home/ubuntu/scratch/tongchen/spCas9/results/{num}/{num}.csv'
    with open(csv_file_path, 'r') as file:
        lines = file.readlines()
    sequence = lines[1].split(',')[1].strip()
    sequence_list.append(sequence)
    length_list.append(len(sequence))
    
assert len(length_list) == len(score_list)

data = {
    "No.": number_list,
    "Length": length_list,
    "TM align score of non-RuvC": score_list,
    "Sequence": sequence_list
}
df = pd.DataFrame(data)

df.to_excel('test_data.xlsx', index=False)