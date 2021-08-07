import click
import trifl

@click.command()
@click.argument('analysis-params-file', type=click.Path(exists=True))
def _analyze(analysis_params_file):
    trifl.analyze(analysis_params_file)
