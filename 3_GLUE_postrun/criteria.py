import os
import json


class Criteria(object):
    def __init__(self, config=None):
        # Define GLUE settings
        self.config = config
        self.id = int(self.config["GLUE"]["criteria_id"])
        self.get_behavioral_criteria()
        self.get_full_behavioral_criteria()

    def get_behavioral_criteria(self):
        """Read the GLUE criteria from JSON file"""
        GLUE_criteria_path = os.path.join(
            self.config["PATHS"]["GLUE_criteria_config"],
            f"GLUE_criteria_{self.id}.json",
        )

        with open(GLUE_criteria_path, "r") as file:
            self.criteria = json.load(file)

    def get_full_behavioral_criteria(self):
        self.full_criteria = dict()
        i = 0

        for key, criterion in self.criteria.items():
            if criterion["metric"] == "season_transition":
                transitions = [
                    "dry2wet_start",
                    "dry2wet_end",
                    "wet2dry_start",
                    "wet2dry_end",
                ]
                for transition in transitions:
                    self.full_criteria.update(
                        {
                            f"{i}": {
                                "metrics_fullname": f"SeasonTrans of Soil {transition}",
                                "threshold": criterion["threshold"],
                                "operation": criterion["operation"],
                            }
                        }
                    )

            else:
                self.full_criteria.update(
                    {
                        f"{i}": {
                            "metrics_fullname": f'{criterion["metric"]} on {criterion["variable_to_analyze"]}',
                            "threshold": criterion["threshold"],
                            "operation": criterion["operation"],
                        }
                    }
                )

            i = +1
