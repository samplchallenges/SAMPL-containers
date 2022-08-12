import json
import os.path
import random
import statistics
import argparse


def score_evaluation(docking_score_prediction, docking_score_answerkey,
                     compound_binds_prediction, compound_binds_answerkey):
    # Actually, load the files, calc RMSE, etc
    result = 1 if compound_binds_prediction == compound_binds_answerkey else 0
    print(f"is_correct {result}")

def score_submissionrun(scores):
    scores_file = scores
    with open(scores_file, "r") as fp:
        scores_dict = json.load(fp)

    is_corrects = [score_dict["is_correct"] for score_dict in scores_dict]

    ct_correct = 0
    for ic in is_corrects:
        if int(ic) == 0:
            ct_correct += 1

    total = len(is_corrects)

    print(f"num_correct {ct_correct}")
    print(f"percent_correct {ct_correct*100/total}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('command')
    parser.add_argument('--docking_score_prediction')
    parser.add_argument('--docking_score_answerkey')
    parser.add_argument('--compound_binds_prediction')
    parser.add_argument('--compound_binds_answerkey')
    parser.add_argument('--scores')
    

    args = parser.parse_args()

    if args.command == 'score-evaluation':
        score_evaluation(args.docking_score_prediction, args.docking_score_answerkey,
                     args.compound_binds_prediction, args.compound_binds_answerkey)

    elif args.command == 'score-submissionrun':
        score_submissionrun(args.scores)

    else:
        print("Invalid command: choose from \'score-evaluation\' or \'score-submissionrun\'.")
