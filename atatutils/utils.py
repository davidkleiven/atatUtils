from typing import List, Tuple, Union

ASE_SELECT_COND = Tuple[str, str, Union[str, int]]

COMPARATORS = ['>=', '<=', '>', '<', '=']

def find_comparator(query: str) -> str:
    for c in COMPARATORS:
        if c in query:
            return c
    raise ValueError("Did not find any of {COMPARATORS} in the query")

def str2select_cond(query: str) -> List[ASE_SELECT_COND]:
    select_conds = []
    for item in query.split(','):
        comparator = find_comparator(query)
        key, value = item.split(comparator)
        
        if value.lower() == 'true':
            value = 1
        elif value.lower() == 'false':
            value = 0
        
        try:
            value = int(value)
        except:
            pass

        select_conds.append((key, comparator, value))
    return select_conds
