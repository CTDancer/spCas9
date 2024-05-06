import random
import os
from ColabFold.colabfold.utils import setup_logging
from ColabFold.colabfold.batch import get_queries, run, set_model_type
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from pathlib import Path
import time

def mutate():
    initial_seq = 'HVQADARVAGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHRLEESFLVEEDKKHERHPIFGNIVDEVAYHEKYPTIYHLRKKLVDSTDKADLRLIYLALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQLVQTYNQLFEENPINASGVDAKAILSARLSKSRRLENLIAQLPGEKKNGLFGNLIALSLGLTPNFKSNFDLAEDAKLQLSKDTYDDDLDNLLAQIGDQYADLFLAAKNLSDAILLSDILRVNTEITKAPLSASMIKRYDEHHQDLTLLKALVRQQLPEKYKEIFFDQSKNGYAGYIDGGASQEEFYKFIKPILEKMDGTEELLVKLNREDLLRKQRTFDNGSIPHQIHLGELHAILRRQEDFYPFLKDNREKIEKILTFRIPYYVGPLARGNSRFAWMTRKSEETITPWNFEEVVDKGASAQSFIERMTNFDKNLPNEKVLPKHSLLYEYFTVYNELTKVKYVTEGMRKPAFLSGEQKKAIVDLLFKTNRKVTVKQLKEDYFKKIECFDSVEISGVEDRFNASLGTYHDLLKIIKDKDFLDNEENEDILEDIVLTLTLFEDREMIEERLKTYAHLFDDKVMKQLKRRRYTGWGRLSRKLINGIRDKQSGKTILDFLKSDGFANRNFMQLIHDDSLTFKEDIQKAQVSGQGDSLIMINLGSPAMKRIVYTVKVVDETVGRKEIVREEWREYQTTQKGQKNSRERMKRIEEGIKELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDHIVPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMKNYWRQLLNAKLITQRKFDNLTKAERGGLSELDKAGFIKRQLVNTRQITKHVDQNLHSMNTNNDYLPEVKVGLKYVMDFRKDYKVEINNYHQAHTEYLNAVVGTAIKKYPKLESEFVYDYKVDVRKMIAKSEQEIGKATAKYFFGLLLPGGLGLEVQTGGFSKESILPKRNSDKLIARKKDWDPKKYGGFDSPTVAYSVLVVAKVEKGKSKKLKSVKELLGITIMERSSFEKNPIDFLEAKGYKEVKKDLIIKLPKYSLFELENGRKRMLASAGELQKGNELALPSKYVNFLYLASHYEKLKGSPEDNEQKQLFVEQHKHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGAPAAFKYFDTTIDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD'
    # intervals = [[1, 25], [689, 739], [892, 958]]
    # non_ruvc_regions = [[26, 688], [740, 891], [1019, 1288]]
    
    intervals = [(1, 9), (673, 713), (866, 977)]
    non_ruvc_regions = [(10, 672), (714, 865), (978, 1247)]
    
    for i in range(209, 1000):
        print(f"########################## Round {i} ##########################")
        start = time.time()
        mutated_seq, mutated_intervals, mutated_non_ruvc_regions = random_mutation(initial_seq, intervals, non_ruvc_regions)
    
        mutated_seq, directory, previous = af2(mutated_seq, i)
        if previous == True:
            mutated_seq = initial_seq
            mutated_intervals = intervals
            mutated_non_ruvc_regions = non_ruvc_regions
        
        pdb_files = find_pdb_files(directory)
        pdb_files = [file for file in pdb_files if 'ruvc' not in file]
        assert len(pdb_files) == 1
        print(f"intervals: {intervals}")
        print(f"mutated_intervals: {mutated_intervals}")
        print(f"non_ruvc_regions: {non_ruvc_regions}")
        print(f"mutated_non_ruvc_regions: {mutated_non_ruvc_regions}")
        tm_score = split(pdb_files[0], mutated_non_ruvc_regions)
        
        if tm_score >= 0.95:
            job_name = directory.split('/')[-1]
            with open("good_results.txt", 'a') as f:
                f.write(f"{job_name} | {tm_score}\n")
            
            initial_seq = mutated_seq
            intervals = mutated_intervals
            non_ruvc_regions = mutated_non_ruvc_regions
        else:
            print("Not better")
            
        print(f"TIME = {time.time() - start}")
        
def find_pdb_files(directory):
    pdb_files = []
    for file in os.listdir(directory):
        if file.endswith(".pdb"):
            pdb_files.append(os.path.join(directory, file))
    return pdb_files


def get_pdb_data_from_file(file):
    s1 = get_structure(file)
    chain1 = next(s1.get_chains())
    coords1, seq1 = get_residue_data(chain1)
    return((coords1, seq1))


def compute_score(pdb1, pdb2):
    coords1, seq1 = get_pdb_data_from_file(pdb1)
    coords2, seq2 = get_pdb_data_from_file(pdb2)

    res = tm_align(coords1, coords2, seq1, seq2)
    
    assert res.tm_norm_chain1 == res.tm_norm_chain2

    print(f"TM-score= {res.tm_norm_chain2}")
    
    return res.tm_norm_chain2


def tm_score(pdb_file, non_ruvc_file, ruvc_file):
    # print("### Total")
    # compute_score(pdb_file, "original.pdb")
    
    # print("\n### Non-RuvC part")
    return compute_score(non_ruvc_file, "original-non_ruvc.pdb")

    # print("\n### RuvC part")
    # compute_score(ruvc_file, "original-ruvc.pdb")
   
        
def split(pdb_file, non_ruvc_regions):
    # print(f"\n-------------{pdb_file.split('/')[-1]}-------------")

    ruvc_file = pdb_file[:-4] + '-ruvc.pdb'
    non_ruvc_file = pdb_file[:-4] + '-non_ruvc.pdb'

    with open(pdb_file, "r") as f:
        lines = f.readlines()

    with open(non_ruvc_file, "w") as f1:
        with open(ruvc_file, "w") as f2:
            for line in lines:
                if 'A' not in line:
                    f1.write(line)
                    f2.write(line)
                else:
                    index = int(line[22:27].strip())
                    if index in range(non_ruvc_regions[0][0], non_ruvc_regions[0][1]+1) or \
                        index in range(non_ruvc_regions[1][0], non_ruvc_regions[1][1]+1) or \
                            index in range(non_ruvc_regions[2][0], non_ruvc_regions[2][1]+1):
                                f1.write(line)
                    else:
                        f2.write(line)
                        
    with open (ruvc_file, "r") as f:
        lines = f.readlines()

    temp = 0
    for i in range(len(lines)-1, 1, -1):
        if 'A' in lines[i]:
            another_line = lines[i][:26]
            middle = str(int(another_line[7:12].strip())+1)
            latter = another_line[17:]
            new_line = "TER    " + middle + "      " + latter + '\n'
            temp = i
            break

    new_lines = lines[:temp+1]
    new_lines.append(new_line)
    for i in range(temp+1, len(lines)):
        new_lines.append(lines[i])

    with open(ruvc_file, "w") as f:
        f.writelines(new_lines)
        
    return tm_score(pdb_file=pdb_file, non_ruvc_file=non_ruvc_file, ruvc_file=ruvc_file)

def af2(query_sequence, i):
    jobname = str(i)
    template_mode = "none"
    num_relax = 0
    
    # def check(folder):
    #     if os.path.exists(folder):
    #         return False
    #     else:
    #         return True
    # if not check(jobname):
    #     n = 0
    # while not check(f"{jobname}_{n}"): n += 1
    
    jobname = f"results/{jobname}"
    if os.path.exists(jobname):
        with open(os.path.join(jobname,f"{jobname.split('/')[-1]}.csv"), "r") as f:
            query_sequence = f.readlines()[-1].split(',')[-1]
        print(f"Previous sequence = {query_sequence}")
        return query_sequence, jobname, True
    
    # make directory to save results
    os.makedirs(jobname, exist_ok=True)
    
    # save queries
    queries_path = os.path.join(jobname, f"{jobname.split('/')[-1]}.csv")
    with open(queries_path, "w") as text_file:
        text_file.write(f"id,sequence\n{jobname},{query_sequence}")
        
    custom_template_path = None
    use_templates = False
    
    print("jobname",jobname)
    print("sequence",query_sequence)
    print("length",len(query_sequence.replace(":","")))
    
    msa_mode = "mmseqs2_uniref_env"
    pair_mode = "unpaired_paired"
    
    a3m_file = os.path.join(jobname,f"{jobname}.a3m")
        
    model_type = "auto"
    num_recycles = "12"
    recycle_early_stop_tolerance = "0.5"
    relax_max_iterations = 200
    pairing_strategy = "greedy"
    max_msa = "auto"
    num_seeds = 1
    use_dropout = False
    num_recycles = None if num_recycles == "auto" else int(num_recycles)
    recycle_early_stop_tolerance = None if recycle_early_stop_tolerance == "auto" else float(recycle_early_stop_tolerance)
    if max_msa == "auto": max_msa = None
    
    save_all = False
    save_recycles = False
    
    try:
        K80_chk = os.popen('nvidia-smi | grep "Tesla K80" | wc -l').read()
    except:
        K80_chk = "0"
    pass
    if "1" in K80_chk:
        print("WARNING: found GPU Tesla K80: limited to total length < 1000")
        if "TF_FORCE_UNIFIED_MEMORY" in os.environ:
            del os.environ["TF_FORCE_UNIFIED_MEMORY"]
        if "XLA_PYTHON_CLIENT_MEM_FRACTION" in os.environ:
            del os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"]
            
    result_dir = jobname
    log_filename = os.path.join(jobname,"log.txt")
    setup_logging(Path(log_filename))
    
    queries, is_complex = get_queries(queries_path)
    model_type = set_model_type(is_complex, model_type)
    
    if "multimer" in model_type and max_msa is not None:
        use_cluster_profile = False
    else:
        use_cluster_profile = True
        
    results = run(
        queries=queries,
        result_dir=result_dir,
        use_templates=use_templates,
        custom_template_path=custom_template_path,
        num_relax=num_relax,
        msa_mode=msa_mode,
        model_type=model_type,
        num_models=1,
        num_recycles=num_recycles,
        relax_max_iterations=relax_max_iterations,
        recycle_early_stop_tolerance=recycle_early_stop_tolerance,
        num_seeds=num_seeds,
        use_dropout=use_dropout,
        model_order=[1],
        is_complex=is_complex,
        data_dir=Path("./ColabFold"),
        keep_existing_results=False,
        rank_by="auto",
        pair_mode=pair_mode,
        pairing_strategy=pairing_strategy,
        stop_at_score=float(100),
        # prediction_callback=prediction_callback,
        dpi=150,
        zip_results=False,
        save_all=save_all,
        max_msa=max_msa,
        use_cluster_profile=use_cluster_profile,
        # input_features_callback=input_features_callback,
        save_recycles=save_recycles,
        user_agent="colabfold/google-colab-main",
    )
    
    return query_sequence, jobname, False
    

def update_intervals(intervals, mutation_position):
    # 更新区间：如果突变位置在区间内，调整该区间和后续区间的位置
    updated_intervals = []
    for start, end in intervals:
        if mutation_position <= end:
            end -= 1
        if mutation_position <= start:
            start -= 1
        updated_intervals.append((start, end))
    return updated_intervals

def select_interval(intervals):
    # Calculate lengths of each interval
    lengths = [b - a for a, b in intervals]

    # Calculate the total length
    total_length = sum(lengths)

    # Calculate the probability for each interval
    probabilities = [length / total_length for length in lengths]

    # Randomly select an interval based on the probabilities
    selected_interval = random.choices(intervals, weights=probabilities, k=1)

    return selected_interval[0]

def random_mutation(sequence, intervals, non_ruvc_regions):
    # 氨基酸单字母代码
    amino_acids = 'ARNDCEQGHILKMFPSTWYV'

    # 随机选择一个区间
    selected_interval = select_interval(intervals)
    start, end = selected_interval

    # 在选定的区间内随机选择一个位置
    mutation_position = random.randint(start, end)
    print(f"Mutation position = {mutation_position}")

    # 决定是替换还是删除
    if random.choice(['replace', 'delete']) == 'replace':
        print("REPLACE")
        # 替换操作：随机选择一个氨基酸进行替换
        replacement_aa = random.choice(amino_acids)
        while replacement_aa == sequence[mutation_position]:
            replacement_aa = random.choice(amino_acids)
        
        sequence = sequence[:mutation_position - 1] + replacement_aa + sequence[mutation_position:]
    else:
        print("DELETE")
        # 删除操作：从序列中删除该位置的氨基酸
        sequence = sequence[:mutation_position - 1] + sequence[mutation_position:]
        # 更新区间
        intervals = update_intervals(intervals, mutation_position)
        non_ruvc_regions = update_intervals(non_ruvc_regions, mutation_position)

    return sequence, intervals, non_ruvc_regions

if __name__ == "__main__":
    mutate()