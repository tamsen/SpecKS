#import git

GLOBAL_git_URL = "https://github.com/tamsen/SpecKS"
GLOBAL_public_version_info = [
    ["Initiated. ", "March 15, 2024", "v1.0.0.0"],
    ["Dev complete, added final shedding module. ", "March 29, 2024", "v1.1.0.0"],
    ["Fixed per_site_evolutionary_distance, to match Tiley's 0.01268182. ", "Apr 3, 2024", "v1.2.0.0"]
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

    def most_recent_comment(self):
        return self.public_comments[-1][0]

    def to_string(self):

        data = ["\n\nGit repo:\t\t" + self.repo_url,
                "Version number:\t" + str(self.version_num),
                "Version notes:\t" + self.most_recent_comment()]

        if self.repo:
            data.append("Git changeset:\t" + self.sha + "\t" + self.when)

        return "\n".join(data)
