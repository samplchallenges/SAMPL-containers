import json
import os.path
import random
import statistics
import argparse


class Classification:
    TRUE_POSITIVE = 1
    TRUE_NEGATIVE = 2
    FALSE_POSITIVE = 3
    FALSE_NEGATIVE = 4

def classify_prediction(answerkey, prediction):
    if answerkey == 1 and prediction == 1:
        return Classification.TRUE_POSITIVE
    elif answerkey == 1 and prediction == 0:
        return Classification.FALSE_NEGATIVE
    elif answerkey == 0 and prediction == 1
        return Classification.FALSE_POSITIVE
    else: #answerkey == 0 and prediction == 0
        return Classification.TRUE_NEGATIVE

def score_evaluation(docking_score_prediction, docking_score_answerkey,
                     compound_binds_prediction, compound_binds_answerkey):
    # Actually, load the files, calc RMSE, etc
    classification = classify_prediction(compound_binds_answerkey, compound_binds_prediction)
    print(f"classification {result}")
    print(f"docking_score {docking_score_prediction}")

def _load_csv(filename):
    result_dict = {}
    with open(filename, "r", encoding="utf8") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            result_dict[row["id"]] = float(row["value"])
    return result_dict

def _write_csv(filename, score_dict):
    with open(filename, "w", encoding="utf8") as fp:
        writer = csv.DictWriter(fp, ["id", "value"])
        writer.writeheader()
        for k, v in score_dict.items():
            writer.writerow({"id": k, "value": v})


def score_batch(output_dir, docking_score_prediction, docking_score_answerkey,
                compound_binds_prediction, compound_binds_answerkey):
    if output_dir is None:
        output_dir = "/mnt/outputs"

    classification_dict = {}
    docking_score_prediction = _load_csv(docking_score_prediction)
    compound_binds_prediction = _load_csv(compound_binds_prediction)
    compound_binds_answerkey = _load_csv(compound_binds_answerkey)

    if len(compound_binds_answerkey) != len(compound_binds_prediction):
        raise ValueError("Predictions and Answers don't match {len(compound_binds_prediction)} vs {len(compound_binds_answerkey)}")
    for key, prediction in compound_binds_prediction.items():
        try:
            answer = compound_binds_answerkey[key]
            classification_dict[key] = classify_prediction(answer, prediction)
        except KeyError:
            raise ValueError("No answer found for key {key}")

    classification_path = os.path.join(output_dir, "classification.csv")
    docking_score_path = os.path.join(output_dir, "docking_score.csv")

    _write_csv(classification_path, classification_dict)
    _write_csv(docking_score_path, docking_score_prediction)

    print(f"classification {classification_path}")
    print(f"docking_score {docking_score_path}")


def score_submissionrun(scores):
    scores_file = scores
    with open(scores_file, "r") as fp:
        scores_dict = json.load(fp)

    classifications = [score_dict["classification"] for score_dict in scores_dict]

    cls_cts = {
        Classification.TRUE_POSITIVE: 0,
        Classification.TRUE_NEGATIVE: 0,
        Classification.FALSE_POSITIVE: 0,
        Classification.FALSE_NEGATIVE: 0,
    }
    for c in classification:
        cls_cts[c] += 1

    print(f"true_positive_ct {cls_cts[Classification.TRUE_POSITIVE]}")
    print(f"true_negative_ct {cls_cts[Classification.TRUE_NEGATIVE]}")
    print(f"false_positive_ct {cls_cts[Classification.FALSE_POSITIVE]}")
    print(f"false_negative_ct {cls_cts[Classification.FALSE_NEGATIVE]}")

    scores = [score_dict["classification"] for score_dict in scores_dict]

    print(f"avg_docking_score {statistics.mean(scores)}")
    print(f"stdev_docking_score {statistics.stdev(scores)}")
    print(f"min_docking_score {min(scores)}")
    print(f"max_docking_score {max(scores)}")


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
