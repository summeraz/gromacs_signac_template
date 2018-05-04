"""Define the project's workflow logic."""
from flow import FlowProject
from flow import staticlabel
import environment  # Custom environment definition


class ProjectAnalysis(FlowProject):

    @staticlabel()
    def density_calculated(job):
        return job.isfile('rho.txt') and job.isfile('rho.pdf')

    def __init__(self, *args, **kwargs):
        super(ProjectAnalysis, self).__init__(*args, **kwargs)

        self.add_operation(
            name='calc_density',
            cmd='python src/analysis-operations.py calc_density {job._id}',
            post=[self.density_calculated])

if __name__ == '__main__':
    ProjectAnalysis().main()
