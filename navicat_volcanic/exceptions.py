#!/usr/bin/env python


class InputError(Exception):
    """Raised when there is an error in the input."""

    pass


class MissingDataError(Exception):
    """Raised when too many missing energies are found in the input. Setting an inputter might fix this."""

    pass


class ChoiceError(InputError):
    """Raised when a choice (a descriptor, a run mode) is needed and none was given while nobody can be asked
    (stdin is not a terminal), or when the given choice is not valid. The command line exits with status 2."""

    pass


class DescriptorChoiceError(ChoiceError):
    """Raised when a descriptor must be chosen but no choice was given and nobody can be asked (stdin is not a
    terminal), or when the given choice is not valid. ``candidates`` lists the candidate descriptors: (column, r2)
    pairs from user_choose_1_dv / user_choose_2_dv, or {"cols", "idx", "r2", "mae"} dicts from choose_1_dv /
    choose_2_dv."""

    def __init__(self, message, candidates=()):
        super().__init__(message)
        self.candidates = list(candidates)
