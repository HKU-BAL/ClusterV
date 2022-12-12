import sys
from importlib import import_module
REPO_NAME="clusterV"
VERSION="1.1"

CV_SCRIPTS_FOLDER="cv"

cv_folder = [
    "ClusterV",
    "run_initial_call",
    "run_ClusterV",
    "get_consensus",
    "update_vcf_af",
    "tag_read",
    "filter_large_indel",
    "get_low_v_db"
]


def directory_for(submodule_name):
    if submodule_name in cv_folder:
        return CV_SCRIPTS_FOLDER
    return ""


def print_help_messages():
    from textwrap import dedent
    print(dedent("""\
        {0} submodule invocator:
            Usage: python cv.py [submodule] [options of the submodule]
            Version {1}
        Available cv submodules:\n{2}
        """.format(
            REPO_NAME,
            VERSION,
            "\n".join("          - %s" % submodule_name for submodule_name in cv_folder),
        )
    ))


def main():
    if len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print_help_messages()
        sys.exit(0)

    if sys.argv[1] == "-v" or sys.argv[1] == "--version":
        print("ClusterV: version {0}".format(VERSION))
        sys.exit(0)
    print("ClusterV: version {0}".format(VERSION))

    submodule_name = sys.argv[1]
    if (
        submodule_name not in cv_folder
    ):
        sys.exit("[ERROR] Submodule %s not found." % (submodule_name))

    directory = directory_for(submodule_name)
    submodule = import_module("%s.%s" % (directory, submodule_name))

    sys.argv = sys.argv[1:]
    sys.argv[0] += (".py")

    # Note: need to make sure every submodule contains main() method
    submodule.main()
    sys.exit(0)

if __name__ == "__main__":
    main()
