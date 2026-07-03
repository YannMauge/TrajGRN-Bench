from importlib.metadata import version, PackageNotFoundError
from flecs import cell_population, data, decay, intervention, mutation, production, sets, trajectory, utils


__all__ = [
    cell_population,
    data,
    decay,
    intervention,
    mutation,
    production,
    sets,
    trajectory,
    utils,
]

try:
    __version__ = version("flecs")
except PackageNotFoundError:
    __version__ = "unknown"