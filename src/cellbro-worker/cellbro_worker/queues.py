from cellbro_db import models

# def queue_paper_summary(paper: models.Paper, priority: int = 6) -> None:
#     from .tasks import create_paper_summary
#     create_paper_summary.apply_async(args=[paper.id], priority=priority)  # type: ignore


# def queue_paper_keywords_extraction(paper: models.Paper, priority: int = 5) -> None:
#     from .tasks import extract_paper_keywords
#     extract_paper_keywords.apply_async(args=[paper.id], priority=priority)  # type: ignore

