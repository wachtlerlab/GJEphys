from setuptools import setup, find_packages
setup(
    name="GJEphys",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=find_packages(exclude=["^\."]),
    exclude_package_data={'': ["Readme.md"]},
    install_requires=["numpy>=1.11.2",
                      "matplotlib>=1.5.3",
                      "scipy>=0.18.1",
                      "pandas>=0.19.0",
                      "seaborn>=0.7.1",
                      # 0.6 had problems with importing Spike2 Files
                      "neo==0.5",
                      "nixio>=0.1.3",
                      "pylatex",
                      "xlrd>=1.0.0",
                      "openpyxl>=2.4.5",
                      "scikit-learn>=0.18.1",
                      "requests>=2.14.2",
                      "traceback2>=1.4"],

    python_requires=">=2.7"

    )