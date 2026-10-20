"""Research-only parser for SNPedia's documented genoset criteria grammar.

No eval, Python execution, or network requests. Unknown calls stay unknown,
including under negation. This module is not a clinical interpretation engine.
"""
import re


class UnsupportedRule(ValueError):
    pass


def parse_rule(text):
    text = re.sub(r'<!--.*?-->', '', text, flags=re.S)
    text = re.sub(r'(?m)^\s*#.*$', '', text).strip()
    if len(text) > 200_000:
        raise UnsupportedRule('Rule too large')
    tokens = re.findall(r'[A-Za-z0-9_-]+|[(),;]|\S', text)
    index = 0

    def take(expected=None):
        nonlocal index
        if index >= len(tokens):
            raise UnsupportedRule('Unexpected end')
        token = tokens[index].lower()
        index += 1
        if expected is not None and token != expected:
            raise UnsupportedRule('Unexpected token')
        return token

    def expression(depth=0):
        if depth > 100:
            raise UnsupportedRule('Rule too deep')
        name = take()
        if re.fullmatch(r'gs\d+', name):
            return ('genoset', name)
        if re.fullmatch(r'(rs|i)\d+', name):
            take('(')
            alleles = [take().upper()]
            if index < len(tokens) and tokens[index] == ';':
                take(';')
                alleles.append(take().upper())
            take(')')
            if any(not re.fullmatch(r'[ACGTID-]+', a) for a in alleles):
                raise UnsupportedRule('Unsupported allele')
            return ('call', name, tuple(sorted(alleles)))
        if name not in {'and', 'or', 'not', 'atleast'}:
            raise UnsupportedRule('Unsupported function or atom')
        take('(')
        threshold = None
        if name == 'atleast':
            value = take()
            if not value.isdigit():
                raise UnsupportedRule('Invalid threshold')
            threshold = int(value)
            take(',')
        children = [expression(depth + 1)]
        while index < len(tokens) and tokens[index] == ',':
            take(',')
            if index < len(tokens) and tokens[index] == ')':
                break
            children.append(expression(depth + 1))
        take(')')
        if threshold is not None and not 0 <= threshold <= len(children):
            raise UnsupportedRule('Threshold outside argument count')
        return (name, threshold, children)

    tree = expression()
    if index != len(tokens):
        raise UnsupportedRule('Trailing content')
    return tree


def evaluate(tree, calls, genoset_lookup):
    """Return True, False or None (unknown); calls use SNPedia orientation."""
    kind = tree[0]
    if kind == 'call':
        observed = calls.get(tree[1])
        if observed is None:
            return None
        desired = tree[2]
        if len(desired) == 1:
            return desired[0] in observed
        return tuple(sorted(observed)) == desired
    if kind == 'genoset':
        return genoset_lookup(tree[1])
    values = [evaluate(child, calls, genoset_lookup) for child in tree[2]]
    yes, unknown = values.count(True), values.count(None)
    if kind == 'and':
        return False if False in values else None if unknown else True
    if kind == 'or':
        return True if yes else None if unknown else False
    if kind == 'not':
        return False if yes else None if unknown else True
    threshold = tree[1]
    return True if yes >= threshold else False if yes + unknown < threshold else None


def evaluate_all(criteria, calls):
    trees, errors, results = {}, {}, {}
    cyclic = set()
    for key, text in criteria.items():
        try:
            trees[key] = parse_rule(text)
        except UnsupportedRule as error:
            errors[key] = str(error)

    def resolve(key, active=frozenset()):
        if key in results:
            return results[key]
        if key in active:
            cyclic.update(active)
            return None
        if key not in trees:
            return None

        def dependency(other):
            # Live rules contain forward references despite the older grammar
            # page's lower-number restriction. Resolve dependencies, not order.
            return resolve(other, active | {key})

        value = evaluate(trees[key], calls, dependency)
        results[key] = None if key in cyclic else value
        return results[key]

    for key in criteria:
        results[key] = resolve(key)
    return results, errors
