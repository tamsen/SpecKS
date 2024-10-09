"""version.py
Handles and reports version information
"""
#import git

GLOBAL_git_URL = "https://github.com/tamsen/SpecKS"
GLOBAL_public_version_info = [
    ["Initiated. ", "March 15, 2024", "v1.0.0.0"],
    ["Dev complete, added final shedding module. ", "March 29, 2024", "v1.1.0.0"],
    ["Fixed per_site_evolutionary_distance, to match Tiley's 0.01268182. ", "Apr 3, 2024", "v1.2.0.0"],
    ["Made gene_div_time_distribution_parameters be a polyploid-level param, not sim-level. ",
     "May 2, 2024", "v1.3.0.0"],
    ["Adding capacity to model segmental allo/autopolyploids. ",
     "May 2, 2024", "v1.4.0.0"],
    ["Retooling with emphasis on 2D continuum. And changed SPC terminology to DIV. ",
    "Aug 6, 2024", "v1.5.0.0"],
    ["Adding reviewer 3 minor suggestions. Affects comments & exposes varaibles, no change to results.",
    "Aug 6, 2024", "v1.5.0.1"]
]



class version_info:

    def __init__(self):

        try:
            #self.repo = git.Repo(search_parent_directories=True)
            self.sha = str(self.repo.head.commit)
            self.when = str(self.repo.head.commit.committed_datetime)

        except:
            self.repo = False
            self.sha = False
            self.when = False

        self.repo_url = GLOBAL_git_URL
        self.app_name = GLOBAL_git_URL.split("/")[-1]

        most_recent_update = GLOBAL_public_version_info[-1]
        self.version_num = most_recent_update[-1]
        self.public_comments = most_recent_update
        self.reference = "'Accurate Inference of the Polyploid Continuum using Forward-time Simulations' "+\
            "by T. Dunn and A. Sethuraman, 2024. Draft currently available at https://www.biorxiv.org/content/10.1101/2024.05.17.594724v2"

    def most_recent_comment(self):
        return self.public_comments[-1][0]

    def to_string(self):

        data = ["\n\nGit repo:\t\t" + self.repo_url,
                "Version number:\t" + str(self.version_num),
                "Reference:\t" + self.reference ]

        if self.repo:
            data.append("Git changeset:\t" + self.sha + "\t" + self.when)

        return "\n".join(data)
